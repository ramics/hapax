"""How much signal a single Illumina read carries, measured rather than asserted.

This exists because the read-level detection direction rests on a claim that is
easy to state and easy to get wrong: that at sub-1x coverage the only evidence
available is what one read (or one pair) carries, and that composition-style
features are too weak at that length to survive a one-in-a-million
false-positive budget. The translated-search arm already measured a clean
negative for homology (results/abundance_sweep/REPORT.md). This script measures
the other candidate features' *ceiling* — what they could possibly discriminate
before any classifier is fitted — on real sequence that is already on disk.

Three sections, each answering a question that would otherwise be hand-waved:

1. ORF-likeness. The usual argument is "random DNA has stops, coding DNA does
   not", which is a statement about kilobases. At the 107 nt read length of
   SRR23973193 a read is ~35 codons, so the per-frame stop-free probability is
   (61/64)^35 and six frames get six attempts. Measured here on real genomes and
   on composition-matched shuffles of the same genomes, because the honest
   control for "is this read coding" is not uniform random sequence but sequence
   with the same base composition. If real non-coding and shuffled sequence are
   both frequently stop-free in some frame, the feature cannot separate them at
   this length whatever the model on top.

2. Frame identifiability. A feature that needs the reading frame (codon usage,
   codon-position bias, the RAMICS-style codon HMM) first needs the frame. Here
   the truth is known: reads are drawn inside annotated HIV-1 CDS intervals, so
   the true frame is the draw offset. Reported as the fraction of reads where
   exactly one of six frames is stop-free, i.e. where the frame is recoverable
   from the read alone, against the fraction where several frames survive and
   the feature has to be averaged over an ambiguous frame.

3. The abundance and budget arithmetic. Derived, not simulated, and kept here so
   that the numbers quoted in the write-up come from a committed script. Two
   parts: how many reads an organism at a given relative abundance contributes
   to a run of a given size, and what per-read false-positive rate a method
   needs for those reads not to be buried. Includes the mate-pair conjunction
   model with an explicit correlation term, because the attractive version of
   the paired-end argument ("require both mates, the false-positive rate
   squares") is only true under independence, and background mates come from the
   same molecule of the same organism, so they are positively correlated. The
   size of that correlation is the quantity that decides whether the pair is
   worth two reads or barely more than one, and it is not measured anywhere in
   the literature surveyed.

Nothing here fits a classifier or claims a detection method works. It bounds
what the features can do, which is the cheap half of the question and the half
that kills candidates before they are built.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import write_manifest  # noqa: E402

STOPS = ("TAA", "TAG", "TGA")
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}


def read_fasta_one(path: Path) -> str:
    seq = []
    with path.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)


def revcomp(s: str) -> str:
    return "".join(COMP.get(c, "N") for c in reversed(s))


def stop_free_frames(seq: str) -> list[int]:
    """Which of the six frames contain no stop codon.

    Frames 0-2 are forward offsets, 3-5 the reverse complement. A read with no
    stop-free frame cannot be coding in any orientation (barring sequencing
    error); a read with several has not told us which frame it is.
    """
    out = []
    rc = revcomp(seq)
    for i, s in enumerate((seq, seq, seq, rc, rc, rc)):
        off = i % 3
        codons = [s[j:j + 3] for j in range(off, len(s) - 2, 3)]
        if not any(c in STOPS for c in codons):
            out.append(i)
    return out


def gc_fraction(seq: str) -> float:
    n = sum(seq.count(b) for b in "ACGT")
    return (seq.count("G") + seq.count("C")) / n if n else 0.0


def iid_stop_free_prob(gc: float, read_len: int) -> dict:
    """Closed form for the same quantity under an iid base model.

    Kept beside the measurement so that a difference between real sequence and
    the model is visible as a difference, not absorbed into a single number. At
    GC 0.5 the per-codon stop probability is the familiar 3/64; higher GC makes
    stops rarer, so high-GC background is *more* often stop-free and the feature
    gets weaker, not stronger, in exactly the samples with high-GC bacteria.
    """
    pg = pc = gc / 2
    pa = pt = (1 - gc) / 2
    p_stop = pt * pa * pa + pt * pa * pg + pt * pg * pa
    n_codons = read_len // 3
    per_frame = (1 - p_stop) ** n_codons
    return {
        "gc": gc,
        "read_len": read_len,
        "p_stop_codon": p_stop,
        "per_frame_stop_free": per_frame,
        # Six frames are not independent (they share bases), so this is an
        # approximation and is labelled as one wherever it is quoted.
        "any_of_six_approx": 1 - (1 - per_frame) ** 6,
    }


def shuffle_seq(seq: str, rng: random.Random) -> str:
    """Mononucleotide shuffle: preserves base composition, destroys everything
    else. The right null for "is stop-free-ness just GC content", which is the
    first objection to any composition feature."""
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def sample_reads(seq: str, read_len: int, n: int, rng: random.Random) -> list[str]:
    reads = []
    hi = len(seq) - read_len
    if hi <= 0:
        return reads
    while len(reads) < n:
        i = rng.randrange(hi)
        r = seq[i:i + read_len]
        if r.count("N") == 0:
            reads.append(r)
    return reads


def section_orf_likeness(genomes: dict[str, Path], read_lens: list[int],
                         n_reads: int, seeds: list[int]) -> list[dict]:
    rows = []
    for name, path in genomes.items():
        seq = read_fasta_one(path)
        gc = gc_fraction(seq)
        for read_len in read_lens:
            for seed in seeds:
                rng = random.Random(seed)
                shuf = shuffle_seq(seq, rng)
                for source, s in (("real", seq), ("shuffled", shuf)):
                    reads = sample_reads(s, read_len, n_reads, random.Random(seed * 7919 + read_len))
                    if not reads:
                        continue
                    frames = [stop_free_frames(r) for r in reads]
                    any_free = sum(1 for f in frames if f)
                    one_free = sum(1 for f in frames if len(f) == 1)
                    rows.append({
                        "genome": name, "source": source, "gc": round(gc, 4),
                        "read_len": read_len, "seed": seed, "n_reads": len(reads),
                        "frac_any_frame_stop_free": any_free / len(reads),
                        "frac_exactly_one_frame": one_free / len(reads),
                        "mean_stop_free_frames": sum(len(f) for f in frames) / len(reads),
                    })
    return rows


def parse_gb_sequence(path: Path) -> str:
    """The ORIGIN block of the same GenBank record the CDS intervals come from.

    Taking the sequence from the record rather than from a fasta of a related
    accession is not fussiness: NC_001802 is HXB2-derived but is 9,181 nt where
    K03455 is 9,719, so K03455's CDS coordinates land in the wrong place on it.
    An earlier pass did exactly that and reported that only ~21% of reads drawn
    "inside a CDS" had a stop-free true frame, and a mean of 0.62 stop-free
    frames per read at 250 nt - below the 1.0 that is guaranteed when a read
    really is inside a coding sequence, which is what gave the bug away.
    """
    out, in_origin = [], False
    for line in path.read_text().splitlines():
        if line.startswith("ORIGIN"):
            in_origin = True
            continue
        if in_origin:
            if line.startswith("//"):
                break
            out.append(re.sub(r"[^acgtnACGTN]", "", line))
    return "".join(out).upper()


def parse_gb_cds(path: Path) -> list[tuple[int, int]]:
    """Contiguous CDS intervals from a GenBank flatfile, 1-based inclusive.

    Deliberately ignores join(...) features: a spliced CDS has no single frame
    for a read spanning the junction, and including them would smuggle an
    unanswerable case into a measurement about frame recovery.
    """
    out = []
    for line in path.read_text().splitlines():
        m = re.match(r"^     CDS             (\d+)\.\.(\d+)$", line)
        if m:
            out.append((int(m.group(1)), int(m.group(2))))
    return out


def section_frame_identifiability(gb: Path, read_lens: list[int],
                                  n_reads: int, seeds: list[int]) -> list[dict]:
    """Reads drawn inside annotated CDS, so the true frame is known exactly.

    Two numbers matter. Whether the true frame is among the stop-free frames
    tells us the feature is not actively wrong. Whether it is the *only*
    stop-free frame tells us the read identifies its own frame, which is what
    any codon-level scoring downstream requires.
    """
    seq = parse_gb_sequence(gb)
    cds = parse_gb_cds(gb)
    rows = []
    for read_len in read_lens:
        for seed in seeds:
            rng = random.Random(seed * 104729 + read_len)
            drawn = 0
            true_in = 0
            true_only = 0
            n_free_total = 0
            while drawn < n_reads:
                start, end = cds[rng.randrange(len(cds))]
                lo, hi = start - 1, end  # to 0-based half-open
                if hi - lo <= read_len:
                    continue
                i = rng.randrange(lo, hi - read_len)
                r = seq[i:i + read_len]
                if r.count("N"):
                    continue
                # Offset of the first whole codon *within the read*, which is
                # what stop_free_frames indexes. Codon boundaries sit at
                # lo, lo+3, ...; a read starting at i therefore has its first
                # boundary (lo - i) mod 3 bases in. Note the sign: (i - lo) % 3
                # is the read's phase within the CDS, not the frame offset, and
                # using it makes the true frame look stop-free in only ~1/3 of
                # reads - which is how this bug was caught, since a CDS has no
                # in-frame stops and the number has to be ~1.
                true_frame = (lo - i) % 3
                free = stop_free_frames(r)
                drawn += 1
                n_free_total += len(free)
                if true_frame in free:
                    true_in += 1
                if free == [true_frame]:
                    true_only += 1
            rows.append({
                "genome": gb.name, "read_len": read_len, "seed": seed,
                "n_reads": drawn,
                "frac_true_frame_stop_free": true_in / drawn,
                "frac_true_frame_unique": true_only / drawn,
                "mean_stop_free_frames": n_free_total / drawn,
            })
    return rows


def section_abundance(run_sizes: list[int], abundances: list[float],
                      genome_len: int, span_nt: int) -> list[dict]:
    """Reads and coverage contributed by an organism at a given abundance.

    `span_nt` is sequence observed per pair, not insert size: a 107 nt mate 1
    plus a 93 nt trimmed mate 2 observes 200 nt of the molecule however long the
    insert is. Coverage is what the assembly arm measured against, so the rungs
    from that arm (nothing below 1x, 1.68x = 4.2% breadth, 5.03x = 64.5%) can be
    read straight off this table.
    """
    rows = []
    for R in run_sizes:
        for a in abundances:
            pairs = a * R
            rows.append({
                "run_pairs": R, "abundance": a,
                "expected_pairs_from_organism": pairs,
                "expected_reads": 2 * pairs,
                "coverage_x": pairs * span_nt / genome_len,
                # Any feature that needs two reads of the same organism needs
                # this to be comfortably above 1, which is the whole reason
                # grouping dies at the bottom of the abundance range.
                "p_at_least_two_pairs": 1 - math.exp(-pairs) * (1 + pairs) if pairs < 700 else 1.0,
            })
    return rows


def section_budget(fp_per_million: list[float], run_sizes: list[int],
                   true_reads: list[int]) -> list[dict]:
    """What a per-read false-positive rate costs at run scale.

    The quantity a surveillance pipeline actually cares about is the ratio of
    false flags to true flags, so this reports both counts. A rate of 1,489
    background reads per million - the lowest non-zero operating point the
    translated-search arm reached - is included so the gap is stated in the same
    units as everything else.
    """
    rows = []
    for R in run_sizes:
        for f in fp_per_million:
            false_flags = R * 2 / 1e6 * f  # reads, not pairs
            for t in true_reads:
                rows.append({
                    "run_pairs": R, "fp_per_million_reads": f,
                    "expected_false_flags": false_flags,
                    "true_reads_available": t,
                    "false_per_true": false_flags / t if t else float("inf"),
                })
    return rows


def section_pair_conjunction(per_read_sens: list[float], per_read_fp: list[float],
                             rho_true: list[float], rho_fp: list[float]) -> list[dict]:
    """Requiring both mates to pass, under stated correlations between mates.

    Model: each mate passes with marginal probability p, and the two mates'
    outcomes are correlated with coefficient rho, so
        P(both) = p^2 + rho * p * (1 - p).
    rho = 0 is independence - the flattering case, where the false-positive rate
    squares - and rho = 1 is total dependence, where the pair is worth one read
    and the conjunction buys nothing.

    The two correlations are kept separate because the asymmetry between them is
    the entire question. A truly novel organism's two mates share an organism
    and a novelty, so rho_true should be high and the conjunction costs little
    sensitivity. Whether rho_fp is also high depends on *why* background reads
    get flagged. If they are flagged because of the organism they came from -
    a divergent but real astrovirus in the sample - then both mates are flagged
    together, rho_fp is high, and requiring both mates removes almost no false
    positives. If they are flagged for read-specific accidents - low complexity,
    adapter remnants, a chimeric molecule, a quality-collapsed tail - then the
    accidents are largely independent between mates, rho_fp is low, and the
    conjunction is worth orders of magnitude.

    So the paired-end conjunction is favourable exactly when background false
    flags are read-level accidents rather than organism-level properties. That
    is a falsifiable statement about one measurable quantity, and it is not
    measured anywhere in the surveyed literature. The measurement it asks for is
    cheap: one translated-search run with mate identity retained, which the
    existing read-level arm discarded by scoring reads as an unpaired pool.
    """
    rows = []
    for s in per_read_sens:
        for f in per_read_fp:
            for rt in rho_true:
                for rf in rho_fp:
                    both_t = s * s + rt * s * (1 - s)
                    both_f = f * f + rf * f * (1 - f)
                    rows.append({
                        "per_read_sensitivity": s,
                        "per_read_fp": f,
                        "rho_true": rt,
                        "rho_fp": rf,
                        "pair_sensitivity": both_t,
                        "pair_fp": both_f,
                        "fp_reduction_factor": f / both_f if both_f else float("inf"),
                        "sens_cost_factor": both_t / s if s else 0.0,
                        # Figure of merit: how much the false-positive rate
                        # falls per unit of sensitivity given up. Above 1 the
                        # trade is favourable; at 1 the conjunction is pointless.
                        "merit": ((f / both_f) / (s / both_t)) if both_f and both_t else float("inf"),
                    })
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--n-seeds", type=int, default=5,
                    help="repo convention: 5 seeds before a number is quoted")
    ap.add_argument("--n-reads", type=int, default=20000)
    ap.add_argument("--read-lens", default="75,107,150,250")
    ap.add_argument("--community", type=Path,
                    default=Path("/mnt/project-files/stress/community/"))
    ap.add_argument("--hiv-gb", type=Path,
                    default=Path("/mnt/project-files/ref/K03455.gb"))
    ap.add_argument("--out", type=Path, default=Path("results/readlevel_features"))
    args = ap.parse_args()

    seeds = list(range(args.seed, args.seed + args.n_seeds))
    read_lens = [int(x) for x in args.read_lens.split(",")]
    args.out.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()

    # A spread of GC and of genome type, so that a feature's behaviour is not
    # read off one organism. Names are accessions; see stress/community.
    wanted = {
        "E_coli_NC_000913.3": "NC_000913.3.fa",
        "S_aureus_NC_007795.1": "NC_007795.1.fa",
        "phage_lambda_NC_001416.1": "NC_001416.1.fa",
        "phiX174_NC_001422.1": "NC_001422.1.fa",
        "SARS2_NC_045512.2": "NC_045512.2.fa",
        "HIV1_NC_001802.1": "NC_001802.1.fa",
        "rabies_NC_001542.1": "NC_001542.1.fa",
        "HBV_NC_003977.2": "NC_003977.2.fa",
    }
    genomes = {}
    for name, fn in wanted.items():
        p = args.community / fn
        if p.exists():
            genomes[name] = p
    if not genomes:
        sys.exit(f"no genomes found under {args.community}")

    timings = {}
    t = time.perf_counter()
    orf = section_orf_likeness(genomes, read_lens, args.n_reads, seeds)
    timings["orf_likeness"] = round(time.perf_counter() - t, 1)

    model = [iid_stop_free_prob(gc, L) for L in read_lens
             for gc in (0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65)]

    frame = []
    if args.hiv_gb.exists():
        # Sequence and CDS coordinates both come from the K03455 record itself;
        # see parse_gb_sequence for why a related accession's fasta will not do.
        t = time.perf_counter()
        frame = section_frame_identifiability(args.hiv_gb, read_lens,
                                              max(2000, args.n_reads // 10), seeds)
        timings["frame_identifiability"] = round(time.perf_counter() - t, 1)

    abundance = section_abundance(
        run_sizes=[int(1e6), int(1e7), int(1e8), int(1e9)],
        abundances=[1e-7, 1e-6, 1e-5, 1e-4, 1e-3],
        genome_len=11932,  # rabies NC_001542, the spike the sweep used
        span_nt=200,       # 107 nt mate 1 + 93 nt trimmed mate 2
    )
    budget = section_budget(
        fp_per_million=[0.1, 1.0, 10.0, 100.0, 1489.0, 93478.0],
        run_sizes=[int(1e6), int(1e8)],
        true_reads=[2, 20, 200],
    )
    pair = section_pair_conjunction(
        per_read_sens=[0.1, 0.3, 0.6],
        per_read_fp=[1e-3, 1e-4, 1e-5],
        rho_true=[0.5, 0.75, 0.9],
        rho_fp=[0.0, 0.1, 0.25, 0.5, 0.75, 0.9],
    )

    payload = {
        "orf_likeness": orf,
        "iid_model": model,
        "frame_identifiability": frame,
        "abundance": abundance,
        "budget": budget,
        "pair_conjunction": pair,
    }
    (args.out / "readlevel_features.json").write_text(json.dumps(payload, indent=2))

    timings["total"] = round(time.perf_counter() - t0, 1)
    write_manifest(args.out, seeds, vars(args) | {
        "community": str(args.community), "hiv_gb": str(args.hiv_gb),
        "out": str(args.out),
    }, timings, extra={"genomes": {k: str(v) for k, v in genomes.items()}})
    print(json.dumps({"sections": {k: len(v) for k, v in payload.items()}, "timings": timings}, indent=2))


if __name__ == "__main__":
    main()
