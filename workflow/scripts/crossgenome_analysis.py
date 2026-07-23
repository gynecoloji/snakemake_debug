#!/usr/bin/env python3
"""Step 09 cross-genome check for ONE sample: of the MULTIMAPPED reads, how many map to
BOTH the host and the spike-in genome?

The combined-genome BAM keeps only each read's best location, and exhaustive re-alignment
(bowtie2 -a) is intractable on repeats, so we cannot read "both genomes" off the BAM. Instead we
re-align the multimapped reads SEPARATELY to the host-only and spike-in-only bowtie2 indexes
(default best-alignment mode -> fast, definitive): a read/fragment maps to BOTH iff it aligns to
both; "codominant" iff the two best scores are within --codominant-margin (an equally-good hit in
each genome = genuinely ambiguous, the reads that distort spike-in normalization).

PAIRED-END aware (auto-detected): if the BAM is paired, whole read PAIRS are re-aligned paired-end
(-X 3000 --no-mixed --no-discordant, matching a typical spike-in pipeline) and classified per
FRAGMENT -- this recovers mate-anchored reads that a single-end re-alignment would drop into
"neither". Single-end BAMs use per-read single-end re-alignment.

By default ALL multimapped reads/fragments are used (--fraction 1.0). Set --fraction < 1 to sample.
Memory: only the spike-in-side score dict is held (the small genome); host alignments are STREAMED.
samtools + bowtie2 + stdlib only."""
import argparse, subprocess, os, re, random, shlex

AS = re.compile(r"\tAS:i:(-?\d+)")
YS = re.compile(r"\tYS:i:(-?\d+)")


def q(s):
    return shlex.quote(str(s))


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def is_paired(align):
    p = subprocess.Popen(["samtools", "view", "-F", "0x900", align],
                         stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    ln = p.stdout.readline()
    p.stdout.close()
    try:
        p.terminate()
    except Exception:
        pass
    return bool(ln) and (int(ln.split("\t", 2)[1]) & 1)


# ---------------- single-end path ----------------
def extract_mm_reads(align, frac, seed, cap, fq_path):
    """Mate-tagged single-end FASTQ of multimapped reads. ALL by default; -s + reservoir if sampling."""
    all_mode = (frac >= 1.0) and (not cap or cap <= 0)
    view = ["samtools", "view", "-F", "0x904", "-e", "[XS]", align]
    if not all_mode:
        view = ["samtools", "view", "-s", f"{seed + min(frac, 0.999999):.6f}",
                "-F", "0x904", "-e", "[XS]", align]
    proc = subprocess.Popen(view, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
                            text=True, bufsize=1 << 20)
    n = 0
    if all_mode or not cap:
        with open(fq_path, "w") as o:
            for ln in proc.stdout:
                f = ln.split("\t", 11)
                if len(f) < 11:
                    continue
                mate = "/1" if (int(f[1]) & 64) else "/2"
                o.write(f"@{f[0]}{mate}\n{f[9]}\n+\n{f[10]}\n")
                n += 1
    else:
        rng = random.Random(seed)
        reservoir = []
        seen = 0
        for ln in proc.stdout:
            f = ln.split("\t", 11)
            if len(f) < 11:
                continue
            mate = "/1" if (int(f[1]) & 64) else "/2"
            rec = f"@{f[0]}{mate}\n{f[9]}\n+\n{f[10]}\n"
            seen += 1
            if len(reservoir) < cap:
                reservoir.append(rec)
            else:
                j = rng.randint(0, seen - 1)
                if j < cap:
                    reservoir[j] = rec
        with open(fq_path, "w") as o:
            o.writelines(reservoir)
        n = len(reservoir)
    proc.stdout.close()
    proc.wait()
    return n


def se_scores(fastq, index, threads):
    """qname -> best AS (single-end). Held in memory -> spike-in (small) genome only."""
    proc = subprocess.Popen(
        ["bowtie2", "-x", index, "-U", fastq, "--no-unal", "-p", str(threads), "--reorder"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
    best = {}
    for ln in proc.stdout:
        if ln[0] == "@":
            continue
        f = ln.split("\t", 12)
        if int(f[1]) & 0x900:
            continue
        m = AS.search(ln)
        s = int(m.group(1)) if m else -10 ** 9
        if f[0] not in best or s > best[f[0]]:
            best[f[0]] = s
    proc.stdout.close()
    proc.wait()
    return best


def se_classify(fastq, index, threads, spk, margin):
    """Stream host single-end alignments past the spike-in dict. -> (n_host, both, codom)."""
    proc = subprocess.Popen(
        ["bowtie2", "-x", index, "-U", fastq, "--no-unal", "-p", str(threads), "--reorder"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
    n_host = both = codom = 0
    for ln in proc.stdout:
        if ln[0] == "@":
            continue
        f = ln.split("\t", 12)
        if int(f[1]) & 0x900:
            continue
        n_host += 1
        d = spk.get(f[0])
        if d is not None:
            both += 1
            m = AS.search(ln)
            h = int(m.group(1)) if m else -10 ** 9
            if abs(h - d) <= margin:
                codom += 1
    proc.stdout.close()
    proc.wait()
    return n_host, both, codom


# ---------------- paired-end path ----------------
def extract_mm_pairs(align, frac, seed, cap, r1, r2, threads):
    """FASTQ of complete read PAIRS for multimapped fragments (any mate with XS). Returns n pairs."""
    names = r1 + ".names"
    samp = "" if (frac >= 1.0) else f"-s {seed + min(frac, 0.999999):.6f} "
    sh(f"samtools view {samp}-F 0x904 -e '[XS]' {q(align)} | cut -f1 | sort -u -S 1G > {q(names)}")
    sh(f"samtools view -N {q(names)} -u -F 0x900 {q(align)} | "
       f"samtools collate -O -u -@ {threads} - {q(r1 + '.ct')} | "
       f"samtools fastq -n -1 {q(r1)} -2 {q(r2)} -0 /dev/null -s /dev/null -@ {threads} - "
       f">/dev/null 2>&1")
    os.remove(names)
    with open(r1) as f:
        n = sum(1 for _ in f) // 4
    return n


def _pair_score(ln):
    a = AS.search(ln)
    y = YS.search(ln)
    return (int(a.group(1)) if a else -10 ** 9) + (int(y.group(1)) if y else 0)


def pe_scores(r1, r2, index, threads):
    """qname -> best concordant-pair score (AS+YS). Held in memory -> spike-in genome only."""
    proc = subprocess.Popen(
        ["bowtie2", "-x", index, "-1", r1, "-2", r2, "-p", str(threads),
         "-X", "3000", "--no-mixed", "--no-discordant", "--no-unal", "--reorder"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
    best = {}
    for ln in proc.stdout:
        if ln[0] == "@":
            continue
        f = ln.split("\t", 12)
        flag = int(f[1])
        if (flag & 0x900) or (flag & 0xC3) != 0x43:  # concordant primary FIRST mate (0x1|0x2|0x40)
            continue
        s = _pair_score(ln)
        if f[0] not in best or s > best[f[0]]:
            best[f[0]] = s
    proc.stdout.close()
    proc.wait()
    return best


def pe_classify(r1, r2, index, threads, spk, margin):
    """Stream host concordant pairs past the spike-in dict. -> (n_host_pairs, both, codom)."""
    proc = subprocess.Popen(
        ["bowtie2", "-x", index, "-1", r1, "-2", r2, "-p", str(threads),
         "-X", "3000", "--no-mixed", "--no-discordant", "--no-unal", "--reorder"],
        stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, bufsize=1 << 20)
    n_host = both = codom = 0
    for ln in proc.stdout:
        if ln[0] == "@":
            continue
        f = ln.split("\t", 12)
        flag = int(f[1])
        if (flag & 0x900) or (flag & 0xC3) != 0x43:  # concordant primary FIRST mate only
            continue
        n_host += 1
        d = spk.get(f[0])
        if d is not None:
            both += 1
            if abs(_pair_score(ln) - d) <= margin:
                codom += 1
    proc.stdout.close()
    proc.wait()
    return n_host, both, codom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--align", required=True)              # combined host+spike-in BAM
    ap.add_argument("--sample", required=True)
    ap.add_argument("--host-index", required=True)
    ap.add_argument("--spikein-index", required=True)
    ap.add_argument("--spikein-prefix", default="spikein_")
    ap.add_argument("--paired", default="auto", choices=["auto", "yes", "no"])
    ap.add_argument("--fraction", type=float, default=1.0)
    ap.add_argument("--seed", type=int, default=13)
    ap.add_argument("--cap", type=int, default=0)
    ap.add_argument("--codominant-margin", type=int, default=5)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)

    paired = is_paired(a.align) if a.paired == "auto" else (a.paired == "yes")
    if paired:
        r1, r2 = a.out + ".R1.fq", a.out + ".R2.fq"
        n = extract_mm_pairs(a.align, a.fraction, a.seed, a.cap, r1, r2, a.threads)
        spk = pe_scores(r1, r2, a.spikein_index, a.threads)
        n_host, both, codom = pe_classify(r1, r2, a.host_index, a.threads, spk, a.codominant_margin)
        for p in (r1, r2):
            try:
                os.remove(p)
            except OSError:
                pass
        unit = "fragments"
    else:
        fq = a.out + ".reads.fq"
        n = extract_mm_reads(a.align, a.fraction, a.seed, a.cap, fq)
        spk = se_scores(fq, a.spikein_index, a.threads)
        n_host, both, codom = se_classify(fq, a.host_index, a.threads, spk, a.codominant_margin)
        try:
            os.remove(fq)
        except OSError:
            pass
        unit = "reads"

    host_only = n_host - both
    spk_only = len(spk) - both
    neither = n - host_only - both - spk_only
    with open(a.out, "w") as o:
        o.write("sample\tn_reads\tboth\tboth_codominant\thost_only\tspike_only\tneither\n")
        o.write(f"{a.sample}\t{n}\t{both}\t{codom}\t{host_only}\t{spk_only}\t{neither}\n")
    pct = lambda x: x / max(n, 1) * 100
    print(f"[09x] {a.sample} ({unit}, {'paired' if paired else 'single'}): tested={n} | "
          f"BOTH genomes={both} ({pct(both):.2f}%, codominant {codom}={pct(codom):.2f}%) "
          f"host_only={host_only} ({pct(host_only):.2f}%) spike_only={spk_only} ({pct(spk_only):.2f}%) "
          f"neither={neither} ({pct(neither):.2f}%)")


if __name__ == "__main__":
    main()
