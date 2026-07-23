#!/usr/bin/env python3
"""Build a single self-contained HTML report from whatever step outputs exist
under --out-dir. One section per step (01–08 + whole-library), theme-aware,
colorblind-safe palette, hover tooltips. Pure stdlib; no external JS/CSS.

    python build_report.py --out-dir <OUT> --out <OUT>/report.html \
        --host "Homo sapiens" --contaminant-label dm6.fa
"""
import argparse, csv, glob, json, os, re
from datetime import datetime, timezone


def rd(p):
    if not os.path.exists(p):
        return []
    with open(p) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fnum(x, d=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


def is_lowcomplex(seq):
    """Crude low-complexity / tandem-repeat flag: few distinct 3-mers for the length."""
    if len(seq) < 12:
        return True
    kmers = {seq[i:i + 3] for i in range(len(seq) - 2)}
    return len(kmers) / (len(seq) - 2) < 0.5


def parse_fasta_counts(path):
    """top_unaligned.fasta headers look like >top07_totalcount=63_len=61 ; return
    {id: {'n':count, 'seq':sequence}} keyed by the leading topNN id."""
    out = {}
    if not os.path.exists(path):
        return out
    cur = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                h = line[1:]
                tid = h.split("_", 1)[0]
                m = re.search(r"totalcount=(\d+)", h)
                cur = tid
                out[tid] = {"n": int(m.group(1)) if m else 0, "seq": ""}
            elif cur is not None:
                out[cur]["seq"] += line.strip()
    return out


def build_data(O, host, contam_label):
    steps = {}

    # ---- 01 alignment rate ----
    r = [x for x in rd(f"{O}/01_alignment_summary/alignment_summary.tsv")
         if x.get("alignment_rate_pct") not in ("", "NA", None)]
    if r:
        rows = sorted(({"sample": x["sample"],
                        "rate": fnum(x["alignment_rate_pct"], 0.0),
                        "total": int(fnum(x.get("total_reads"), 0))} for x in r),
                      key=lambda d: d["rate"])
        hi = max(d["rate"] for d in rows)
        lo = min(d["rate"] for d in rows)
        steps["01"] = {"rows": rows, "lo": lo, "hi": hi,
                       "low": [d["sample"] for d in rows if d["rate"] < hi - 10]}

    # ---- 02 signature composition ----
    r = rd(f"{O}/02_sequence_signatures/signature_fractions.tsv")
    if r:
        cats = ["bact_16S", "illumina_adapter", "polyA", "host_repeat", "other"]
        rows = []
        for x in r:
            row = {"sample": x["sample"]}
            for c in cats:
                row[c] = round(fnum(x.get(c), 0.0) * 100, 3)
            rows.append(row)
        # order by "other" ascending so the informative fractions read left-to-right
        rows.sort(key=lambda d: d["other"])
        steps["02"] = {"rows": rows, "cats": cats}

    # ---- 03 BLAST identity of top sequences ----
    counts = parse_fasta_counts(f"{O}/03_blast_top_sequences/top_unaligned.fasta")
    if counts:
        hits = {}
        for x in rd(f"{O}/03_blast_top_sequences/blast_best_hits.tsv"):
            q = x.get("query", "")
            tid = q.split("_", 1)[0]
            hits.setdefault(tid, x.get("subject", ""))
        seqs = []
        for tid, info in counts.items():
            subj = hits.get(tid, "")
            low = subj.lower()
            if subj:
                if "mitochond" in low or "d-loop" in low:
                    cat, label = "human_mito", "Human mtDNA D-loop"
                elif "drosophila" in low:
                    cat, label = "fly", "Drosophila melanogaster"
                elif "homo sapiens" in low or "human" in low:
                    cat, label = "human", "Homo sapiens"
                else:
                    cat, label = "other_hit", subj[:48]
                call = subj[:70]
            else:
                if is_lowcomplex(info["seq"]):
                    cat, label = "repeat", "simple / satellite repeat"
                    call = "no nt hit — low-complexity / tandem repeat"
                else:
                    cat, label = "nohit", "no BLAST hit"
                    call = "no nt hit returned"
            seqs.append({"id": tid, "n": info["n"], "cat": cat,
                         "label": label, "call": call, "seq": info["seq"]})
        seqs.sort(key=lambda d: -d["n"])
        steps["03"] = {"seqs": seqs, "blast_run": bool(hits)}

    # ---- 04 random-read BLAST species ----
    sp = rd(f"{O}/04_blast_random/random_blast_species.tsv")
    if sp:
        steps["04"] = {"species": [{"sp": x["species"], "reads": int(fnum(x["reads"], 0)),
                                     "pct": fnum(x.get("pct_of_hits"), 0.0)} for x in sp]}

    # ---- 05 contaminant genome screen ----
    r = rd(f"{O}/05_contaminant_genome_screen/contaminant_alignment.tsv")
    if r:
        rows = sorted(({"sample": x["sample"],
                        "pct": fnum(x["pct_reads_mapping_contaminant"], 0.0)} for x in r),
                      key=lambda d: -d["pct"])
        steps["05"] = {"rows": rows, "label": contam_label or "contaminant genome"}

    # ---- 06 Kraken2 ----
    r = rd(f"{O}/06_kraken2/kraken_summary.tsv")
    if r:
        rows = []
        for x in r:
            total = fnum(x.get("total_reads"), 0.0) or 1.0
            classified = fnum(x.get("pct_classified"), 0.0)
            m = re.search(r"Homo sapiens\((\d+)\)", x.get("top_species", ""))
            human = round(int(m.group(1)) / total * 100, 3) if m else 0.0
            human = min(human, classified)
            rows.append({"sample": x["sample"],
                         "uncls": round(100 - classified, 3),
                         "human": human,
                         "other": round(max(classified - human, 0.0), 3),
                         "bacteria": fnum(x.get("pct_bacteria"), 0.0)})
        rows.sort(key=lambda d: -d["uncls"])
        species = [{"sp": x["species"], "n": int(fnum(x["total_reads_across_samples"], 0))}
                   for x in rd(f"{O}/06_kraken2/top_species_overall.tsv")]
        steps["06"] = {"rows": rows, "species": species[:10]}

    # ---- 07 custom-sequence mapping ----
    cm = rd(f"{O}/07_custom_sequences/custom_mapping.tsv")
    if cm:
        ps = rd(f"{O}/07_custom_sequences/custom_per_sequence.tsv")
        seqcols = [c for c in (ps[0].keys() if ps else []) if c != "sample"]
        clean = lambda n: re.sub(r"_5'.*$", "", n)
        rows = [{"sample": x["sample"], "pct": fnum(x["pct_mapping_custom"], 0.0),
                 "mapped": int(fnum(x.get("reads_mapping_custom"), 0)),
                 "total": int(fnum(x.get("total_unaligned_reads"), 0))} for x in cm]
        perseq = {x["sample"]: {c: int(fnum(x.get(c), 0)) for c in seqcols} for x in ps}
        steps["07"] = {"rows": rows, "seqs_raw": seqcols,
                       "seqs": [clean(c) for c in seqcols], "perseq": perseq}

    # ---- 08 top-organism (auto-fetched RefSeq) genome screen ----
    to = rd(f"{O}/08_top_organism/top_organism_alignment.tsv")
    if to and to[0].get("organism") and to[0].get("organism") != "NONE":
        rows = [{"sample": x["sample"], "pct": fnum(x["pct_reads_mapping_top_organism"], 0.0)} for x in to]
        steps["08"] = {"rows": sorted(rows, key=lambda r: -r["pct"]),
                       "organism": to[0]["organism"], "accession": to[0].get("accession", "")}

    # ---- 09 multimapped reads ----
    mm = rd(f"{O}/09_multimapped_reads/multimap_summary.tsv")
    if mm:
        rows = [{"sample": x["sample"], "method": x.get("method", ""),
                 "pct_mm": fnum(x["pct_multimapped"], 0.0),
                 "pct_host": fnum(x.get("pct_host"), 0.0),
                 "pct_spikein": fnum(x.get("pct_spikein"), 0.0),
                 "pct_host_mm": fnum(x.get("pct_host_multimapped"), 0.0),
                 "pct_spikein_mm": fnum(x.get("pct_spikein_multimapped"), 0.0)} for x in mm]
        spike = any(r["pct_spikein"] > 0 for r in rows)
        loci = [{"contig": x["contig"], "reads": int(fnum(x.get("multimapped_reads"), 0)),
                 "pct": fnum(x.get("pct_of_multimapped"), 0.0)}
                for x in rd(f"{O}/09_multimapped_reads/top_multimap_loci.tsv")]
        xg = [{"sample": x["sample"], "pct_both": fnum(x.get("pct_both"), 0.0),
               "pct_codom": fnum(x.get("pct_both_codominant"), 0.0),
               "pct_host_only": fnum(x.get("pct_host_only"), 0.0),
               "pct_spike_only": fnum(x.get("pct_spike_only"), 0.0),
               "pct_neither": fnum(x.get("pct_neither"), 0.0),
               "n": int(fnum(x.get("n_reads"), 0))}
              for x in rd(f"{O}/09_multimapped_reads/crossgenome_summary.tsv")]
        steps["09"] = {"rows": sorted(rows, key=lambda r: -r["pct_mm"]), "spike": spike,
                       "loci": loci[:12], "method": rows[0]["method"] if rows else "",
                       "xg": sorted(xg, key=lambda r: -r["pct_both"])}

    # ---- whole-library composition: host / top contamination / others ---------
    # Use the precise bowtie2 genome alignment for the dominant contaminant — the
    # auto top-organism screen (step 08) if present, else the step-05 contaminant
    # genome — and lump everything else in the unaligned pool as "others".
    if "01" in steps:
        aln = {d["sample"]: d["rate"] for d in steps["01"]["rows"]}
        top, label = {}, ""
        if "08" in steps:
            top = {d["sample"]: d["pct"] for d in steps["08"]["rows"]}
            label = steps["08"]["organism"]
        elif "05" in steps:
            top = {d["sample"]: d["pct"] for d in steps["05"]["rows"]}
            label = contam_label or "contaminant genome"
        lib = []
        for s, host in aln.items():
            un = max(100 - host, 0.0)
            tf = min(top.get(s, 0.0) / 100, 1.0)
            lib.append({"sample": s, "host": round(host, 2),
                        "contam": round(un * tf, 2), "other": round(un * (1 - tf), 2)})
        lib.sort(key=lambda r: r["host"])        # worst-aligning (most contaminated) first
        if lib:
            steps["lib"] = {"rows": lib, "label": label}

    return {"host": host or "host",
            "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
            "steps": steps}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--host", default="host")
    ap.add_argument("--contaminant-label", default="contaminant genome")
    a = ap.parse_args()

    data = build_data(a.out_dir, a.host, a.contaminant_label)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    html = HEAD + "\n<script>\nconst DATA = " + json.dumps(data) + ";\n" + JS + "\n</script>\n"
    with open(a.out, "w") as f:
        f.write(html)
    print(f"[report] wrote {a.out} ({len(data['steps'])} step sections)")


# ------------------------------------------------------------------ CSS + shell
HEAD = r"""<title>Unaligned-read diagnostics — report</title>
<style>
  :root{
    --bg:#f5f5f2;--panel:#fff;--panel-2:#fafaf8;--ink:#1b2027;--ink-2:#565c64;--muted:#8b9099;
    --line:#e5e4de;--grid:#edece7;--accent:#3a4a63;
    --c-human:#0072B2;--c-fly:#D55E00;--c-sat:#CC79A7;--c-noise:#009E73;--c-uncls:#c3c5c1;
    --shadow:0 1px 2px rgba(20,24,30,.05),0 6px 20px rgba(20,24,30,.05);
  }
  @media (prefers-color-scheme:dark){:root{
    --bg:#111213;--panel:#1c1d1c;--panel-2:#191a19;--ink:#e9e9e5;--ink-2:#b1b3ad;--muted:#7e817b;
    --line:#2b2c29;--grid:#242522;--accent:#9fb2d0;
    --c-human:#2E86C0;--c-fly:#D96D28;--c-sat:#C173A8;--c-noise:#1F9A78;--c-uncls:#4b4e4a;
    --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 26px rgba(0,0,0,.35);}}
  :root[data-theme="dark"]{
    --bg:#111213;--panel:#1c1d1c;--panel-2:#191a19;--ink:#e9e9e5;--ink-2:#b1b3ad;--muted:#7e817b;
    --line:#2b2c29;--grid:#242522;--accent:#9fb2d0;
    --c-human:#2E86C0;--c-fly:#D96D28;--c-sat:#C173A8;--c-noise:#1F9A78;--c-uncls:#4b4e4a;
    --shadow:0 1px 2px rgba(0,0,0,.4),0 8px 26px rgba(0,0,0,.35);}
  :root[data-theme="light"]{
    --bg:#f5f5f2;--panel:#fff;--panel-2:#fafaf8;--ink:#1b2027;--ink-2:#565c64;--muted:#8b9099;
    --line:#e5e4de;--grid:#edece7;--accent:#3a4a63;
    --c-human:#0072B2;--c-fly:#D55E00;--c-sat:#CC79A7;--c-noise:#009E73;--c-uncls:#c3c5c1;
    --shadow:0 1px 2px rgba(20,24,30,.05),0 6px 20px rgba(20,24,30,.05);}
  *{box-sizing:border-box}
  body{background:var(--bg);color:var(--ink);margin:0;font-size:16px;line-height:1.55;
    font-family:ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
    -webkit-font-smoothing:antialiased}
  .mono{font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace}
  .wrap{max-width:1000px;margin:0 auto;padding:clamp(20px,4vw,52px) clamp(16px,4vw,40px) 80px}
  .eyebrow{font-size:12px;letter-spacing:.14em;text-transform:uppercase;color:var(--accent);font-weight:650;margin:0 0 14px}
  h1{font-family:Georgia,"Iowan Old Style",serif;font-weight:600;font-size:clamp(28px,5vw,46px);
    line-height:1.05;letter-spacing:-.015em;margin:0 0 16px;text-wrap:balance}
  .lede{font-size:clamp(15px,1.8vw,18px);color:var(--ink-2);max-width:66ch;margin:0 0 26px}
  .tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:14px;margin:0 0 6px}
  .tile{background:var(--panel);border:1px solid var(--line);border-radius:14px;padding:17px 18px 15px;
    box-shadow:var(--shadow);position:relative;overflow:hidden}
  .tile::before{content:"";position:absolute;left:0;top:0;bottom:0;width:3px;background:var(--bar,var(--accent))}
  .tile .k{font-size:12.5px;color:var(--muted);margin:0 0 6px}
  .tile .v{font-size:clamp(24px,3.6vw,32px);font-weight:680;letter-spacing:-.02em;line-height:1}
  .tile .s{font-size:12.5px;color:var(--ink-2);margin-top:7px}
  section{margin-top:44px}
  .sec-head{display:flex;align-items:baseline;gap:12px;margin:0 0 4px;flex-wrap:wrap}
  .sec-tag{font-family:ui-monospace,monospace;font-size:12px;font-weight:600;color:var(--panel);
    background:var(--accent);padding:2px 8px;border-radius:6px}
  h2{font-size:clamp(19px,2.4vw,24px);font-weight:640;letter-spacing:-.01em;margin:0;text-wrap:balance}
  .sec-sub{color:var(--ink-2);font-size:14.5px;max-width:74ch;margin:10px 0 18px}
  .card{background:var(--panel);border:1px solid var(--line);border-radius:16px;box-shadow:var(--shadow);
    padding:clamp(16px,2.2vw,24px);margin-top:16px}
  .card h3{font-size:14px;font-weight:640;margin:0 0 2px}
  .card .cap{font-size:12.5px;color:var(--muted);margin:0 0 16px}
  .rows{display:flex;flex-direction:column;gap:10px}
  .row{display:grid;grid-template-columns:150px 1fr;align-items:center;gap:14px}
  .row .lab{font-size:13px;color:var(--ink-2);text-align:right;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  .track{position:relative;background:var(--grid);border-radius:6px;height:26px;display:flex}
  .fill{height:100%;border-radius:5px;min-width:2px;display:flex;align-items:center;justify-content:flex-end;transition:filter .12s}
  .fill:hover{filter:brightness(1.08)}
  .fill .bl{font-size:12px;font-weight:640;color:#fff;padding:0 9px;white-space:nowrap;text-shadow:0 1px 1px rgba(0,0,0,.28)}
  .val-out{position:absolute;left:calc(var(--w) + 8px);top:50%;transform:translateY(-50%);font-size:12px;font-weight:640;color:var(--ink-2);white-space:nowrap}
  .seg{height:100%;display:flex;align-items:center;justify-content:center;overflow:hidden}
  .seg:first-child{border-radius:5px 0 0 5px}.seg:last-child{border-radius:0 5px 5px 0}
  .seg+.seg{margin-left:2px}
  .seg .sl{font-size:11px;font-weight:640;color:#fff;text-shadow:0 1px 1px rgba(0,0,0,.3);padding:0 4px;white-space:nowrap}
  .axis{display:grid;grid-template-columns:150px 1fr;gap:14px;margin-top:12px}
  .axis .ticks{position:relative;height:16px;border-top:1px solid var(--line)}
  .axis .t{position:absolute;top:0;font-size:11px;color:var(--muted);transform:translateX(-50%)}
  .axis .t::before{content:"";position:absolute;top:-1px;left:50%;height:4px;width:1px;background:var(--line)}
  .legend{display:flex;flex-wrap:wrap;gap:15px;margin:2px 0 18px}
  .lg{display:flex;align-items:center;gap:7px;font-size:12.5px;color:var(--ink-2)}
  .sw{width:12px;height:12px;border-radius:3px;flex:none}
  .tbl-wrap{overflow-x:auto;margin-top:18px;border:1px solid var(--line);border-radius:12px}
  table{border-collapse:collapse;width:100%;font-size:13px;min-width:480px}
  th,td{text-align:left;padding:9px 14px;border-bottom:1px solid var(--line);white-space:nowrap}
  thead th{font-size:11px;letter-spacing:.05em;text-transform:uppercase;color:var(--muted);font-weight:600;background:var(--panel-2)}
  tbody tr:last-child td{border-bottom:none}
  td.n{text-align:right;font-variant-numeric:tabular-nums}
  .dot{display:inline-block;width:8px;height:8px;border-radius:50%;margin-right:7px;vertical-align:middle}
  .callout{display:flex;gap:14px;align-items:flex-start;background:var(--panel-2);border:1px solid var(--line);
    border-left:3px solid var(--c-fly);border-radius:12px;padding:15px 18px;margin-top:18px}
  .callout p{margin:0;font-size:13.5px;color:var(--ink-2)}.callout b{color:var(--ink)}
  .note{font-size:12.5px;color:var(--muted);margin-top:14px;max-width:78ch}
  footer{margin-top:52px;padding-top:20px;border-top:1px solid var(--line);color:var(--muted);font-size:12.5px}
  footer code{font-family:ui-monospace,monospace;color:var(--ink-2);font-size:11.5px}
  #tip{position:fixed;z-index:50;pointer-events:none;opacity:0;transition:opacity .1s;background:var(--ink);
    color:var(--bg);font-size:12px;padding:7px 10px;border-radius:8px;box-shadow:0 6px 24px rgba(0,0,0,.28);max-width:280px;line-height:1.4}
  #tip .tt{font-weight:700;margin-bottom:2px}
  @media (max-width:560px){.row,.axis{grid-template-columns:96px 1fr}}
  @media (prefers-reduced-motion:reduce){*{transition:none!important}}
</style>
<div class="wrap" id="root"></div>
<div id="tip"></div>"""


# ------------------------------------------------------------------ JS renderer
JS = r"""
const root = document.getElementById('root'), tip = document.getElementById('tip');
const S = DATA.steps;
const C = { human:'var(--c-human)', fly:'var(--c-fly)', sat:'var(--c-sat)', noise:'var(--c-noise)', uncls:'var(--c-uncls)' };
const short = s => s.replace(/^GSF\d+-/, '');
const fmtN = n => n.toLocaleString();
function el(tag, cls, html){ const e=document.createElement(tag); if(cls)e.className=cls; if(html!=null)e.innerHTML=html; return e; }
function bindTip(e, html){
  e.addEventListener('mousemove', ev=>{ tip.innerHTML=html; tip.style.opacity='1';
    let x=ev.clientX+14,y=ev.clientY+14; const r=tip.getBoundingClientRect();
    if(x+r.width>innerWidth-8)x=ev.clientX-r.width-14; if(y+r.height>innerHeight-8)y=ev.clientY-r.height-14;
    tip.style.left=x+'px'; tip.style.top=y+'px'; });
  e.addEventListener('mouseleave', ()=>tip.style.opacity='0');
}
function ticks(max, step, suf){
  const ax=el('div','axis'), sp=el('div'), tk=el('div','ticks');
  for(let v=0; v<=max+1e-6; v+=step){ const t=el('div','t', v+(suf||'')); t.style.left=(v/max*100)+'%'; tk.appendChild(t); }
  ax.appendChild(sp); ax.appendChild(tk); return ax;
}
// single-series horizontal bars
function bars(data, {max, color, label, value, outfmt, tip:tipFn}){
  const box=el('div','rows');
  data.forEach(d=>{
    const w=Math.max(value(d)/max*100, 0);
    const row=el('div','row'); row.appendChild(el('div','lab', label(d)));
    const track=el('div','track'), fill=el('div','fill');
    fill.style.width=w+'%'; fill.style.background=(typeof color==='function')?color(d):color;
    const txt=outfmt?outfmt(d):value(d);
    if(w>16){ fill.appendChild(el('span','bl', txt)); }
    else { const o=el('span','val-out', txt); o.style.setProperty('--w', w+'%'); track.appendChild(o); }
    track.appendChild(fill); if(tipFn) bindTip(fill, tipFn(d));
    row.appendChild(track); box.appendChild(row);
  });
  return box;
}
// horizontal bars with a nested overlay (a subset drawn on top, from the same left origin)
function overlayBars(data, {max, baseColor, overColor, value, over, label, outfmt, tip:tipFn}){
  const box=el('div','rows');
  data.forEach(d=>{
    const wb=Math.max(value(d)/max*100,0), wo=Math.max(over(d)/max*100,0);
    const row=el('div','row'); row.appendChild(el('div','lab', label(d)));
    const track=el('div','track');
    const base=el('div','fill'); base.style.width=wb+'%'; base.style.background=baseColor;
    const txt=outfmt?outfmt(d):value(d);
    if(wb>16){ base.appendChild(el('span','bl', txt)); }
    else { const o=el('span','val-out', txt); o.style.setProperty('--w', wb+'%'); track.appendChild(o); }
    const ov=el('div','fill'); ov.style.position='absolute'; ov.style.left='0'; ov.style.top='0';
    ov.style.width=wo+'%'; ov.style.background=overColor; ov.style.borderRadius='5px';
    track.appendChild(base); track.appendChild(ov);
    if(tipFn){ bindTip(base, tipFn(d)); bindTip(ov, tipFn(d)); }
    row.appendChild(track); box.appendChild(row);
  });
  return box;
}
// stacked-to-100 horizontal bars
function stacked(data, segs, tipFn){
  const box=el('div','rows');
  data.forEach(d=>{
    const row=el('div','row'); row.appendChild(el('div','lab', short(d.sample)));
    const track=el('div','track');
    segs.forEach(([key,col,name])=>{
      const v=d[key]; if(!(v>0)) return;
      const seg=el('div','seg'); seg.style.width=v+'%'; seg.style.background=col;
      if(v>=8) seg.appendChild(el('span','sl', v.toFixed(v<10?1:0)+'%'));
      bindTip(seg, tipFn(d, name, v)); track.appendChild(seg);
    });
    row.appendChild(track); box.appendChild(row);
  });
  return box;
}
function legend(items){
  const lg=el('div','legend');
  items.forEach(([col,txt])=>{ const i=el('span','lg'); const s=el('span','sw'); s.style.background=col;
    i.appendChild(s); i.appendChild(document.createTextNode(txt)); lg.appendChild(i); });
  return lg;
}
function card(title, cap){ const c=el('div','card'); c.appendChild(el('h3',null,title)); if(cap)c.appendChild(el('p','cap',cap)); return c; }
function section(tag, title, sub){
  const s=el('section'); const h=el('div','sec-head');
  h.appendChild(el('span','sec-tag', tag)); h.appendChild(el('h2', null, title));
  s.appendChild(h); if(sub) s.appendChild(el('p','sec-sub', sub)); return s;
}
function table(cols, rows){
  const w=el('div','tbl-wrap'), t=el('table');
  const th=el('thead'), htr=el('tr');
  cols.forEach(c=>{ const e=el('th', c.n?'n':null, c.t); htr.appendChild(e); }); th.appendChild(htr); t.appendChild(th);
  const tb=el('tbody');
  rows.forEach(r=>{ const tr=el('tr'); r.forEach((cell,i)=>{ const td=el('td', cols[i].n?'n':null); td.innerHTML=cell; tr.appendChild(td); }); tb.appendChild(tr); });
  t.appendChild(tb); w.appendChild(t); return w;
}

// ---------------- header ----------------
root.appendChild(el('p','eyebrow', 'Low-alignment / contamination diagnostics · ' + DATA.host));
root.appendChild(el('h1', null, 'Where the unaligned reads come from'));
root.appendChild(el('p','lede',
  'Interactive per-step report generated by the diagnostics workflow. Each section below corresponds to one step of the pipeline; hover any bar for exact values.'));
const tiles=el('div','tiles');
function tile(k,v,s,bar){ const t=el('div','tile'); if(bar)t.style.setProperty('--bar',bar);
  t.appendChild(el('p','k',k)); t.appendChild(el('div','v',v)); if(s)t.appendChild(el('p','s',s)); return t; }
if(S['01']) tiles.appendChild(tile('Alignment rate', S['01'].lo.toFixed(1)+'–'+S['01'].hi.toFixed(1)+'%',
  S['01'].rows.length+' samples', C.human));
if(S['05']){ const p=S['05'].rows.map(r=>r.pct); tiles.appendChild(tile('Unaligned → '+S['05'].label,
  Math.min(...p).toFixed(0)+'–'+Math.max(...p).toFixed(0)+'%', 'of unaligned reads', C.fly)); }
if(S['06']){ const u=S['06'].rows.map(r=>r.uncls); tiles.appendChild(tile('Kraken2 unclassified',
  Math.min(...u).toFixed(0)+'–'+Math.max(...u).toFixed(0)+'%', 'of unaligned reads', C.uncls)); }
if(tiles.children.length) root.appendChild(tiles);

// ---------------- 01 alignment ----------------
if(S['01']){
  const s=section('01','Alignment rate per sample',
    'Overall mapped fraction to the host reference (samtools flagstat). Bars below the best by >10% are flagged.');
  const c=card('% reads aligned to host', 'Hover for read counts.');
  const hi=S['01'].hi;
  c.appendChild(bars(S['01'].rows, { max:100, color:d=> d.rate<hi-10 ? C.fly : C.human,
    label:d=>short(d.sample), value:d=>d.rate, outfmt:d=>d.rate.toFixed(1)+'%',
    tip:d=>'<div class="tt">'+short(d.sample)+'</div>'+d.rate.toFixed(2)+'% aligned<br>'+fmtN(d.total)+' reads' }));
  c.appendChild(ticks(100,25,'%')); s.appendChild(c);
  s.appendChild(table([{t:'Sample'},{t:'Aligned %',n:1},{t:'Reads',n:1}],
    S['01'].rows.map(d=>[short(d.sample), d.rate.toFixed(2), fmtN(d.total)])));
  root.appendChild(s);
}

// ---------------- 02 signatures ----------------
if(S['02']){
  const cats=[['bact_16S',C.noise,'Bacterial 16S'],['illumina_adapter',C.fly,'Illumina adapter'],
    ['polyA',C.sat,'poly-A/T'],['host_repeat',C.human,'Host repeat (Alu)'],['other','var(--c-uncls)','Other / unclassified']];
  const s=section('02','Motif signatures of unaligned reads',
    'Fraction of sampled unaligned reads matching diagnostic motifs. A large "other" means the cause is not in this motif set (see BLAST / Kraken2 / contaminant screen).');
  const c=card('Composition of unaligned reads (%)');
  c.appendChild(legend(cats.map(x=>[x[1],x[2]])));
  c.appendChild(stacked(S['02'].rows, cats,
    (d,name,v)=>'<div class="tt">'+short(d.sample)+' · '+name+'</div><b>'+v.toFixed(2)+'%</b> of unaligned reads'));
  s.appendChild(c); root.appendChild(s);
}

// ---------------- 03 BLAST ----------------
if(S['03']){
  const CAT={ human_mito:[C.human,'Human mtDNA D-loop'], human:[C.human,'Homo sapiens'],
    fly:[C.fly,'Drosophila'], other_hit:[C.noise,'Other BLAST hit'],
    repeat:[C.sat,'Simple / satellite repeat — no hit'], nohit:['var(--c-uncls)','No BLAST hit'] };
  const seqs=S['03'].seqs, max=Math.max(...seqs.map(d=>d.n));
  const present=[...new Set(seqs.map(d=>d.cat))];
  const s=section('03','Identity of the most-duplicated unaligned sequences',
    'Top duplicated unaligned sequences pooled across samples' + (S['03'].blast_run?', BLASTed against NCBI nt.':' (BLAST not run — showing abundance and low-complexity flag).'));
  const c=card('Abundance of each top sequence, by identity', 'Bar = total duplicate count across samples. Hover for the call.');
  c.appendChild(legend(present.map(k=>[CAT[k][0],CAT[k][1]])));
  c.appendChild(bars(seqs, { max, color:d=>CAT[d.cat][0], label:d=>d.id, value:d=>d.n, outfmt:d=>fmtN(d.n),
    tip:d=>'<div class="tt">'+d.id+' · '+fmtN(d.n)+' copies</div>'+d.call }));
  c.appendChild(ticks(Math.ceil(max/300)*300, Math.ceil(max/300)*300/4, ''));
  s.appendChild(c);
  s.appendChild(table([{t:'Seq'},{t:'Count',n:1},{t:'Identity'}],
    seqs.map(d=>['<span class="mono">'+d.id+'</span>', fmtN(d.n),
      '<span class="dot" style="background:'+CAT[d.cat][0]+'"></span>'+d.call])));
  root.appendChild(s);
}

// ---------------- 04 random-read BLAST ----------------
if(S['04'] && S['04'].species.length){
  const sp=S['04'].species, smax=Math.max(...sp.map(d=>d.reads));
  const s=section('04','Random-read BLAST — species discovered',
    'A random (unbiased) sample of unaligned reads BLASTed against NCBI nt. Not biased toward the most-duplicated sequences, so it surfaces organisms that Kraken2 (DB-limited) and the motif screen both miss.');
  const c=card('Species hit by random unaligned reads', 'Bar = number of random reads whose best nt hit is this species. Hover for exact counts.');
  c.appendChild(bars(sp, { max:smax, color:d=> /homo sapiens/i.test(d.sp)?C.human:(/drosophila/i.test(d.sp)?C.fly:C.noise),
    label:d=>d.sp, value:d=>d.reads, outfmt:d=>d.reads+' ('+d.pct.toFixed(0)+'%)',
    tip:d=>'<div class="tt">'+d.sp+'</div><b>'+d.reads+'</b> random reads · '+d.pct.toFixed(1)+'% of hits' }));
  s.appendChild(c);
  s.appendChild(table([{t:'Species'},{t:'Reads',n:1},{t:'% of hits',n:1}],
    sp.map(d=>['<span style="font-style:italic">'+d.sp+'</span>', d.reads, d.pct.toFixed(1)])));
  root.appendChild(s);
}

// ---------------- 05 contaminant ----------------
if(S['05']){
  const p=S['05'].rows.map(r=>r.pct), max=Math.max(100, Math.ceil(Math.max(...p)/25)*25);
  const s=section('05','Unaligned reads realigned to the '+S['05'].label,
    'All of each sample’s unaligned reads realigned to the suspect genome (bowtie2 --very-sensitive-local). A high, sample-varying rate that tracks the alignment deficit indicates real contamination.');
  const c=card('% of unaligned reads mapping to '+S['05'].label);
  c.appendChild(bars(S['05'].rows, { max:100, color:C.fly, label:d=>short(d.sample), value:d=>d.pct,
    outfmt:d=>d.pct.toFixed(1)+'%',
    tip:d=>'<div class="tt">'+short(d.sample)+'</div><b>'+d.pct.toFixed(2)+'%</b> of unaligned reads map to '+S['05'].label }));
  c.appendChild(ticks(100,25,'%')); s.appendChild(c); root.appendChild(s);
}

// ---------------- 06 Kraken2 ----------------
if(S['06']){
  const segs=[['uncls',C.uncls,'Unclassified'],['human',C.human,'Human'],['other',C.noise,'Other (bacterial/viral)']];
  const s=section('06','Kraken2 taxonomic classification',
    'Composition of unaligned reads by Kraken2. A large unclassified fraction means the reads are not represented in the database (e.g. a genome absent from the DB).');
  const c=card('Classification of unaligned reads (stacked to 100%)', 'Hover any segment for exact %.');
  c.appendChild(legend(segs.map(x=>[x[1],x[2]])));
  c.appendChild(stacked(S['06'].rows, segs,
    (d,name,v)=>'<div class="tt">'+short(d.sample)+' · '+name+'</div><b>'+v.toFixed(2)+'%</b> of unaligned reads'));
  s.appendChild(c);
  if(S['06'].species && S['06'].species.length){
    const sp=S['06'].species, smax=Math.sqrt(sp[0].n);
    const c2=card('Top classified species (pooled across samples)', '√-scaled axis so low-count noise is visible beside the dominant organism.');
    c2.appendChild(bars(sp, { max:smax, color:d=> /homo sapiens/i.test(d.sp)?C.human:C.noise,
      label:d=>d.sp, value:d=>Math.sqrt(d.n), outfmt:d=>fmtN(d.n),
      tip:d=>'<div class="tt">'+d.sp+'</div><b>'+fmtN(d.n)+'</b> reads classified' }));
    s.appendChild(c2);
  }
  root.appendChild(s);
}

// ---------------- 07 custom-sequence mapping ----------------
if(S['07'] && S['07'].rows.length){
  const st=S['07'], pmax=Math.max(...st.rows.map(r=>r.pct), 5);
  const PAL=[C.sat, C.fly, C.noise, C.human];
  const s=section('07','Custom-sequence mapping',
    'Unaligned reads mapped (bowtie2) to your custom FASTA — e.g. a transgene / vector construct. Shows how much of the unaligned pool is the construct, and which sequence within it the reads hit.');
  const c=card('% of unaligned reads mapping to the custom sequences', 'Hover for read counts.');
  c.appendChild(bars(st.rows, { max:pmax, color:C.sat, label:d=>short(d.sample), value:d=>d.pct,
    outfmt:d=>d.pct.toFixed(2)+'%',
    tip:d=>'<div class="tt">'+short(d.sample)+'</div><b>'+d.pct.toFixed(2)+'%</b> of unaligned reads<br>'+fmtN(d.mapped)+' / '+fmtN(d.total)+' reads' }));
  s.appendChild(c);
  if(st.seqs.length){
    const c2=card('Which custom sequence the reads hit', 'Composition of the construct-mapped reads per sample (stacked to 100%).');
    c2.appendChild(legend(st.seqs.map((nm,i)=>[PAL[Math.min(i,3)], nm])));
    const stackRows = st.rows.map(r=>{
      const counts=st.seqs_raw.map(nm=> (st.perseq[r.sample]||{})[nm]||0);
      const tot=counts.reduce((a,b)=>a+b,0)||1; const o={sample:r.sample};
      counts.forEach((c,i)=> o['s'+i]=+(c/tot*100).toFixed(2)); return o;
    });
    const segs=st.seqs.map((nm,i)=>['s'+i, PAL[Math.min(i,3)], nm]);
    c2.appendChild(stacked(stackRows, segs,
      (d,name,v)=>'<div class="tt">'+short(d.sample)+' · '+name+'</div><b>'+v.toFixed(1)+'%</b> of construct-mapped reads'));
    s.appendChild(c2);
    s.appendChild(table([{t:'Sample'},{t:'% mapping',n:1},{t:'Reads',n:1}].concat(st.seqs.map(nm=>({t:nm,n:1}))),
      st.rows.map(r=>[short(r.sample), r.pct.toFixed(2), fmtN(r.mapped)]
        .concat(st.seqs_raw.map(nm=>fmtN((st.perseq[r.sample]||{})[nm]||0))))));
  }
  root.appendChild(s);
}

// ---------------- 08 top-organism (auto) genome screen ----------------
if(S['08'] && S['08'].rows.length){
  const st=S['08'], pmax=Math.max(...st.rows.map(r=>r.pct), 5);
  const s=section('08','Top-organism genome screen — '+st.organism,
    'The top non-host organism from the analysis (Kraken2 + random-BLAST) was auto-selected, its RefSeq genome ('+st.accession+') fetched from NCBI, and each sample’s unaligned reads aligned to it — a precise genome-level quantification of the contaminant.');
  const c=card('% of unaligned reads mapping to '+st.organism+' ('+st.accession+')', 'Hover for exact values.');
  c.appendChild(bars(st.rows, { max:Math.min(100,Math.max(pmax,25)), color:C.noise, label:d=>short(d.sample), value:d=>d.pct,
    outfmt:d=>d.pct.toFixed(1)+'%',
    tip:d=>'<div class="tt">'+short(d.sample)+'</div><b>'+d.pct.toFixed(2)+'%</b> of unaligned reads map to '+st.organism }));
  s.appendChild(c); root.appendChild(s);
}

// ---------------- 09 multimapped reads ----------------
if(S['09'] && S['09'].rows.length){
  const st=S['09'];
  const s=section('09','Multimapped reads'+(st.spike?' — spike-in composition':''),
    'Reads that aligned to multiple locations, from the raw aligner output (detection: '+st.method+'). '+
    (st.spike ? 'Reads are split into host vs spike-in genome; the multimapped set is then broken down by genome and by the contigs it lands on — the multimapper distribution.'
              : 'Where they concentrate (below) tells you why — rRNA, repeats, or multi-copy gene families.'));
  if(st.spike){
    const segs=[['pct_host',C.human,'Host genome'],['pct_spikein',C.fly,'Spike-in genome']];
    const c=card('Host vs spike-in genome split (the normalization signal)', 'samtools idxstats over all reads; stacked to 100% of mapped reads. Hover for exact %.');
    c.appendChild(legend(segs.map(x=>[x[1],x[2]])));
    c.appendChild(stacked(st.rows, segs,
      (d,name,v)=>'<div class="tt">'+short(d.sample)+' · '+name+'</div><b>'+v.toFixed(2)+'%</b> of mapped reads'));
    s.appendChild(c);
    const pmax=Math.max(...st.rows.map(r=>r.pct_mm), 2);
    const c2=card('% of mapped reads that are multimapped', 'Aligner-aware ('+st.method+'); for bowtie2 (XS) this reproduces its own ">1 time" alignment rate.');
    c2.appendChild(bars(st.rows, { max:pmax, color:C.sat, label:d=>short(d.sample), value:d=>d.pct_mm,
      outfmt:d=>d.pct_mm.toFixed(2)+'%',
      tip:d=>'<div class="tt">'+short(d.sample)+'</div><b>'+d.pct_mm.toFixed(2)+'%</b> multimapped ('+st.method+')' }));
    s.appendChild(c2);
    // multimapper distribution BY GENOME: what fraction of each genome's own reads is ambiguous
    const byg=[]; st.rows.forEach(d=>{ byg.push({lab:short(d.sample)+' · host',v:d.pct_host_mm,g:'host'});
                                       byg.push({lab:short(d.sample)+' · spike-in',v:d.pct_spikein_mm,g:'spike-in'}); });
    const gmax=Math.max(...byg.map(d=>d.v), 2);
    const c3=card('Multimapping rate within each genome', 'Of the reads assigned to a genome, how many are multimapped. A high spike-in rate means ambiguous spike-in reads that can distort normalization.');
    c3.appendChild(legend([[C.human,'host'],[C.fly,'spike-in']]));
    c3.appendChild(bars(byg, { max:gmax, color:d=>d.g==='host'?C.human:C.fly, label:d=>d.lab, value:d=>d.v,
      outfmt:d=>d.v.toFixed(2)+'%',
      tip:d=>'<div class="tt">'+d.lab+'</div><b>'+d.v.toFixed(2)+'%</b> of '+d.g+' reads multimapped' }));
    s.appendChild(c3);
    s.appendChild(table([{t:'Sample'},{t:'% host',n:1},{t:'% spike-in',n:1},{t:'% multimapped',n:1},{t:'host mm%',n:1},{t:'spike mm%',n:1}],
      st.rows.map(d=>[short(d.sample), d.pct_host.toFixed(2), d.pct_spikein.toFixed(2), d.pct_mm.toFixed(2), d.pct_host_mm.toFixed(2), d.pct_spikein_mm.toFixed(2)])));
    if(st.loci.length){
      const lmax=Math.max(...st.loci.map(d=>d.reads), 1);
      const c4=card('Where the multimappers land (top contigs)', 'Pooled across samples — host chr* and spike-in contigs the ambiguous reads map to.');
      c4.appendChild(bars(st.loci, { max:lmax, color:C.noise, label:d=>d.contig, value:d=>d.reads,
        outfmt:d=>fmtN(d.reads),
        tip:d=>'<div class="tt">'+d.contig+'</div><b>'+fmtN(d.reads)+'</b> multimapped reads · '+d.pct.toFixed(1)+'%' }));
      s.appendChild(c4);
    }
    if(st.xg && st.xg.length){
      const xmax=Math.max(...st.xg.map(d=>d.pct_both), 0.5);
      const c5=card('Multimapped reads that map to BOTH genomes', 'Every multimapped read (paired-end: fragment) re-aligned separately to the host-only and spike-in-only genomes. Bar = % mapping to both; the darker inset = "codominant" (equally-good in each genome → genuinely ambiguous, the reads that distort normalization).');
      c5.appendChild(legend([[C.sat,'maps to both'],[C.fly,'codominant (ambiguous)']]));
      c5.appendChild(overlayBars(st.xg, { max:xmax, baseColor:C.sat, overColor:C.fly,
        value:d=>d.pct_both, over:d=>d.pct_codom, label:d=>short(d.sample),
        outfmt:d=>d.pct_both.toFixed(2)+'%',
        tip:d=>'<div class="tt">'+short(d.sample)+' · n='+fmtN(d.n)+'</div><b>'+d.pct_both.toFixed(2)+'%</b> map to both genomes ('+d.pct_codom.toFixed(2)+'% codominant)' }));
      s.appendChild(c5);
      const xsegs=[['pct_host_only',C.human,'host-only'],['pct_spike_only',C.fly,'spike-only'],
                   ['pct_both',C.sat,'both genomes'],['pct_neither',C.uncls,'neither']];
      const c6=card('Re-alignment outcome of multimapped reads', 'Every multimapped read (paired-end: fragment) re-aligned separately to the host-only and spike-in-only genomes, classified by where it maps. Stacked to 100%.');
      c6.appendChild(legend(xsegs.map(x=>[x[1],x[2]])));
      c6.appendChild(stacked(st.xg, xsegs,
        (d,name,v)=>'<div class="tt">'+short(d.sample)+' · '+name+'</div><b>'+v.toFixed(2)+'%</b> of multimapped reads'));
      s.appendChild(c6);
      s.appendChild(table([{t:'Sample'},{t:'reads/pairs',n:1},{t:'% host-only',n:1},{t:'% spike-only',n:1},{t:'% both',n:1},{t:'% codominant',n:1}],
        st.xg.map(d=>[short(d.sample), fmtN(d.n), d.pct_host_only.toFixed(2), d.pct_spike_only.toFixed(2), d.pct_both.toFixed(3), d.pct_codom.toFixed(3)])));
    }
  } else {
    const pmax=Math.max(...st.rows.map(r=>r.pct_mm), 5);
    const c=card('% of mapped reads that are multimapped', 'Hover for exact values.');
    c.appendChild(bars(st.rows, { max:pmax, color:C.sat, label:d=>short(d.sample), value:d=>d.pct_mm,
      outfmt:d=>d.pct_mm.toFixed(2)+'%',
      tip:d=>'<div class="tt">'+short(d.sample)+'</div><b>'+d.pct_mm.toFixed(2)+'%</b> multimapped ('+st.method+')' }));
    s.appendChild(c);
    if(st.loci.length){
      const lmax=Math.max(...st.loci.map(d=>d.reads), 1);
      const c2=card('Where the multimappers concentrate (top contigs)', 'Pooled across samples — the loci driving multimapping.');
      c2.appendChild(bars(st.loci, { max:lmax, color:C.noise, label:d=>d.contig, value:d=>d.reads,
        outfmt:d=>fmtN(d.reads),
        tip:d=>'<div class="tt">'+d.contig+'</div><b>'+fmtN(d.reads)+'</b> multimapped reads · '+d.pct.toFixed(1)+'%' }));
      s.appendChild(c2);
    }
  }
  root.appendChild(s);
}

// ---------------- whole-library composition ----------------
if(S['lib']){
  const st=S['lib'];
  const clab=(st.label||'contaminant').replace(/\.(fa|fasta|fna)$/i,'');
  const segs=[['host',C.human,'Host (aligned)'],['contam',C.fly,'Top contamination — '+clab],['other','var(--c-uncls)','Others']];
  const s=section('∑','Whole-library composition',
    'Every read split into host (aligned), the top contamination (bowtie2 genome alignment — '+clab+'), and others.');
  const c=card('Share of the library (%)');
  c.appendChild(legend(segs.map(x=>[x[1],x[2]])));
  c.appendChild(stacked(st.rows, segs,
    (d,name,v)=>'<div class="tt">'+short(d.sample)+' · '+name+'</div><b>'+v.toFixed(2)+'%</b> of all reads'));
  s.appendChild(c); root.appendChild(s);
}

// ---------------- footer ----------------
const f=el('footer');
f.appendChild(el('p', null, 'Generated ' + DATA.generated + ' by <code>build_report.py</code> from the step outputs under the diagnostics <code>output_dir</code>. Steps with no output are omitted.'));
root.appendChild(f);
"""

if __name__ == "__main__":
    main()
