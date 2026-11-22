"""Measure, rather than assert, how often two more local copies of the dinucleotide shuffle
fall back to returning their input.

Both carry the Altschul-Erikson connectivity condition that the broken copies lacked, so the
question is not "is it the same bug" but "how often does the retry loop exhaust", which is
the rate their committed nulls would contain unshuffled input at. The rate depends on the
input, so it is measured on sequence of the two lengths those scripts take - reads and
contigs - and reported as that rather than as a general figure.
"""
import random, sys, importlib.util, collections
from pathlib import Path

def load(p, name):
    spec = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m

def seqs_from_fasta(p, limit, minlen=0):
    out, name, buf = {}, None, []
    with open(p) as fh:
        for line in fh:
            if line.startswith(">"):
                if name and len("".join(buf)) >= minlen:
                    out[name] = "".join(buf).upper()
                    if len(out) >= limit: return out
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    return out

def reads_from_fastq(p, limit):
    import gzip
    out = {}
    with gzip.open(p, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                out[str(i)] = line.strip().upper()
                if len(out) >= limit: return out
    return out

S = Path(sys.argv[1])
rng = random.Random(1)
cases = [
    ("residue_structure.dinuc_shuffle", "experiments/residue_structure.py",
     reads_from_fastq(S / "reads/SRR23421544.fastq.gz", 5000), "5,000 reads"),
    ("residue_contig_annotation.dinuc_shuffle", "experiments/residue_contig_annotation.py",
     seqs_from_fasta(S / "asm/contigs_renamed_identified_norrna.fa", 5000, 300),
     "5,000 contigs >= 300 nt"),
]
for label, path, data, what in cases:
    m = load(path, label.split(".")[0] + "_probe")
    f = m.dinuc_shuffle
    identical = badlen = baddi = 0
    for s in data.values():
        try:
            r = f(s, rng)
        except Exception as e:
            print(f"  raised {type(e).__name__}"); continue
        if r == s: identical += 1
        elif len(r) != len(s): badlen += 1
        elif collections.Counter(zip(r, r[1:])) != collections.Counter(zip(s, s[1:])): baddi += 1
    n = len(data)
    print(f"{label} on {what}: returned input unchanged {identical}/{n} "
          f"({100*identical/n:.2f}%), length violations {badlen}, dinucleotide violations {baddi}")

# MEASURED 2026-11-22 on the holdout library, both at seed 1:
#   residue_structure.dinuc_shuffle,          5,000 reads:            0/5000 (0.00%)
#   residue_contig_annotation.dinuc_shuffle,  5,000 contigs >= 300 nt: 0/5000 (0.00%)
# zero length violations and zero dinucleotide-count violations in both.
#
# THIS READING OF ZERO IS CORRECTED AND THE CORRECTION IS THE POINT OF KEEPING IT HERE.
# The `Make hapax fail loudly` thread ran the same function over 7,390 read-length windows,
# with ONE RNG SEED PER CALL the way the script itself calls it rather than one stream shared
# across calls as below, and got 2 returned unchanged - 0.027% - on ordinary 93 and 107 nt
# reads, reproducing a 0.00-0.10% reading from an earlier five-seed run. Their setup is the
# faithful one. My 0 of 5,000 and their 2 of 7,390 are both consistent with a rate near
# 0.03%, which NEITHER sample size pins, so the figure to quote is 0.027% and not zero: see
# the rule that 0 of N bounds a rate at about 1 in N and never establishes zero.
# residue_contig_annotation reads 0 of 7,390 there too, so it agrees with this run.
#
# What the two numbers together say. Both of these copies DO carry the Altschul-Erikson
# connectivity condition that the two genuinely broken copies lacked, so the condition is
# NECESSARY AND NOT SUFFICIENT: a few sequences fail structurally and then fail on every one
# of the 50 retries, which is why more retries would not help. The defect in the broken copies
# was the missing condition; the bare `return s` only made the failure SILENT. Greping for
# the string therefore gives two false positives, and the fix for these two is to raise and
# report the rate rather than to re-derive the walk.
#
# Not re-run for this, by agreement with that thread: 0.027% moves nothing published and the
# direction is conservative (an unshuffled read in the null inflates the null, which weakens
# the test rather than manufacturing a result). The raise-instead-of-fallback change rides
# along the next time each script is touched for its own reason. The measured rates are what
# make deferring it a decision rather than a shrug.
