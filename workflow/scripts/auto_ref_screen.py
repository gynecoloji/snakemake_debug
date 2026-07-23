#!/usr/bin/env python3
"""Step 09 (login node / internet): auto-pick the top non-host organism from the
analysis (Kraken2 + random-BLAST), fetch its RefSeq reference genome from NCBI,
and align each sample's unaligned reads to it -> per-sample % mapping.

Self-contained (fetch + index + align loop) so the dynamically-chosen reference
does not need to be known at DAG-build time."""
import argparse, io, json, os, re, subprocess, urllib.parse, urllib.request, zipfile

from pick_top_organism import pick_organism

API = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
UA = {"User-Agent": "snakemake_debug/1.0"}


def _get(url, timeout=120):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout)


def fetch_reference(name, out_fa):
    """Download a RefSeq genome FASTA for the taxon `name`. Returns accession or None."""
    q = urllib.parse.quote(name)
    acc = None
    for extra in ("&filters.reference_only=true", ""):
        url = f"{API}/genome/taxon/{q}/dataset_report?filters.assembly_source=RefSeq&page_size=5{extra}"
        try:
            reps = json.load(_get(url, 90)).get("reports") or []
        except Exception:
            reps = []
        if reps:
            acc = reps[0].get("accession")
            break
    if not acc:
        return None
    try:
        raw = _get(f"{API}/genome/accession/{acc}/download?include_annotation_type=GENOME_FASTA", 600).read()
        z = zipfile.ZipFile(io.BytesIO(raw))
        fnas = [n for n in z.namelist() if n.endswith((".fna", ".fa", ".fasta"))]
        if not fnas:
            return None
        with open(out_fa, "wb") as o:
            for n in fnas:
                o.write(z.read(n))
        return acc
    except Exception:
        return None


def align_pct(prefix, fq, threads):
    """% of reads in fq that map to the reference (bowtie2 overall alignment rate)."""
    p = subprocess.run(["bowtie2", "-x", prefix, "-U", fq, "-p", str(threads), "--very-sensitive-local"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
    m = re.search(r"([\d.]+)% overall alignment rate", p.stderr)
    return float(m.group(1)) if m else 0.0


def sample_name(path):
    b = os.path.basename(path)
    for suffix in (".all_unaligned.fq.gz", ".unaligned.fq.gz"):
        if b.endswith(suffix):
            return b[:-len(suffix)]
    return os.path.splitext(b)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--host", default="host")
    ap.add_argument("--min-fraction", type=float, default=0.05)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-pick", required=True)
    ap.add_argument("--fastqs", nargs="+", required=True)
    a = ap.parse_args()

    O = a.out_dir
    refdir = f"{O}/08_top_organism/ref"
    os.makedirs(refdir, exist_ok=True)
    for p in (a.out_tsv, a.out_pick):
        os.makedirs(os.path.dirname(p), exist_ok=True)

    name, score, source = pick_organism(f"{O}/06_kraken2/top_species_overall.tsv",
                                        f"{O}/04_blast_random/random_blast_species.tsv",
                                        a.host, a.min_fraction)

    def finish(pick_line, rows, organism="", accession=""):
        with open(a.out_pick, "w") as o:
            o.write(pick_line + "\n")
        with open(a.out_tsv, "w") as o:
            o.write("sample\tpct_reads_mapping_top_organism\torganism\taccession\n")
            for s, pct in rows:
                o.write(f"{s}\t{pct:.4f}\t{organism}\t{accession}\n")

    if not name:
        print("[08] no credible non-host organism passed the threshold; skipping.")
        finish(f"organism\tNONE\nreason\tno non-host organism >= {a.min_fraction} support", [])
        return

    print(f"[08] top non-host organism: {name} (support {score:.2f}, from {source})")
    ref_fa = f"{refdir}/top_organism.fa"
    acc = fetch_reference(name, ref_fa)
    if not acc:
        print(f"[08] could not fetch a RefSeq genome for '{name}'.")
        finish(f"organism\t{name}\nsupport\t{score:.3f}\nsource\t{source}\naccession\tNOT_FOUND", [])
        return
    print(f"[08] fetched {name} genome {acc}; building index + aligning {len(a.fastqs)} samples")

    prefix = f"{refdir}/top_organism"
    subprocess.run(["bowtie2-build", "--threads", str(a.threads), ref_fa, prefix],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    rows = []
    for fq in sorted(a.fastqs):
        rows.append((sample_name(fq), align_pct(prefix, fq, a.threads)))
    rows.sort(key=lambda r: -r[1])

    finish(f"organism\t{name}\nsupport\t{score:.3f}\nsource\t{source}\naccession\t{acc}",
           rows, organism=name, accession=acc)
    print(f"[08] wrote {a.out_tsv} ({len(rows)} samples vs {name} {acc})")


if __name__ == "__main__":
    main()
