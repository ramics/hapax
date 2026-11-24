#!/usr/bin/env python3
"""Build the end-to-end fixture: a tiny library with all three output classes in it.

WHY THIS EXISTS. Nothing in this repo has ever run reads to report. CI installs no
external tool and no index, so every test stops at a unit or a mocked stage, and the
consequence was measured: the first production pass over an untouched library finished
with zero failed stages and built a page that was mostly a census of absent files
(PLAN.md 1.1, 1.3, 1.5). A fixture is what lets the whole command be a test rather than
a thing someone runs by hand on an instance and reads afterwards.

WHAT IT HAS TO CONTAIN, and why it is three genomes and not one. The report has three
output classes on one distance continuum (SPEC 2), so a fixture with only known material
in it would render two thirds of the page as an absence and the test would pass on a tool
that cannot do most of its job. So the library carries, by construction:

  * material that IS identifiable and must be subtracted - MS2, in the toy nucleotide
    index, and E. coli, in the toy classifier library and the toy SILVA;
  * material that is divergent but related - a mutated PMMoV whose protein identity to
    the real PMMoV in the toy protein index is MEASURED here and written into the
    manifest, not assumed, because the mutation procedure decides it and a stated
    divergence nobody checked is how this project has been wrong before;
  * material that is in no index at all, at low depth - the unplaceable arm, which is
    where the residue report and the depth-1 pile come from.

WHAT IT IS NOT. It is not a benchmark and no number from it belongs in `results/`. The
genomes are three, the library is ten thousand pairs, and the indexes are toy indexes
built from the fixture's own sequence: recall against them says nothing about recall
against RefSeq. It answers one question only, which is whether every section of the
report is filled from a file the run wrote. That is why the divergence is recorded with
how it was established rather than used as a result: validation data belongs in tests
and never in the report (SPEC 13.4).

THE ERROR MODEL is the one measured on the discriminating dataset reference-free from
self-overlapping pairs: 0.00068 per base along the body and 0.0019 at the final base
(results/berlin_quality/, CLAUDE.md). A fixture with no errors would let an assembler
look better than it is, and one with an invented error rate would be a number without
provenance.

Usage:
    python tests/fixtures/make_fixture.py --seed 1          # rebuild from NCBI
    python tests/fixtures/make_fixture.py --seed 1 --refs DIR   # rebuild from local copies

The committed fixture is the output of `--seed 1`. Re-running it with the same seed
reproduces the files byte for byte; `manifest.json` carries the seed, the accessions,
the download date, this script's content hash and the measured divergences.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
import subprocess
import sys
import time
import urllib.request
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

#: The three real accessions, and why each one is in the fixture.
#:
#: Real sequence rather than random strings, because a toy index built from random
#: sequence cannot be wrong in any of the ways a real one is: no codon structure, no
#: conserved domains, no repeats, and an aligner's significance statistics are calibrated
#: on real composition. Small real genomes are the cheapest way to keep that.
ACCESSIONS = {
    "ms2": {
        "accession": "NC_001417.2",
        "role": "known virus, in the toy nucleotide index, must be subtracted as tier 1",
        "why": "Enterobacteria phage MS2, 3,569 nt, four non-overlapping ORFs whose roles "
               "are maturation, coat, lysis and replicase. The repo's gene-level "
               "asymmetry work is stated in exactly those roles (SPEC 2.2), so a fixture "
               "carrying them exercises the per-gene ranking rather than only a "
               "whole-contig identity.",
    },
    "pmmov": {
        "accession": "NC_003630.1",
        "role": "divergent-but-related virus: mutated here, the real one in the toy "
                "protein index only",
        "why": "Pepper mild mottle virus, 6,357 nt. A plant virus that is one of the most "
               "abundant RNA viruses in human wastewater, so a wastewater fixture "
               "carrying it is not an arbitrary choice.",
    },
    "ecoli": {
        "accession": "NC_000913.3",
        "role": "bacterial background, in the toy classifier library and the toy SILVA",
        "why": "E. coli K-12 MG1655. A 16 kb slice spanning the rrnB operon, so the "
               "fixture has ribosomal sequence in it: rRNA is about a third of the "
               "discriminating library and the stage that anchors against SILVA has "
               "nothing to do without any.",
    },
}

#: The slice of E. coli that is taken, and why these coordinates.
#:
#: 1-based inclusive, as the GenBank record numbers it. rrnB (16S rrsB, 23S rrlB, 5S
#: rrfB) sits at about 4,166,000-4,172,000, so this window carries a whole rRNA operon
#: plus flanking protein-coding sequence. The whole 4.64 Mb genome would dominate the
#: fixture's size for no gain: the classifier arm needs sequence it can name, not a
#: genome.
ECOLI_SLICE = (4_164_000, 4_180_000)

#: Per-genome coverage in the simulated library, and why these depths.
#:
#: They are not equal, because an equal-abundance library is the one case where nothing
#: in the cascade or the residue report has to deal with the thing it exists for. MS2 and
#: PMMoV are deep enough to assemble on their own; the unplaceable arm sits at 3x so the
#: run has a depth-1 pile, which is where the measured boundary of this tool is
#: (CLAUDE.md: the unplaceable residue sits at median depth 1).
COVERAGE = {"ms2": 300, "pmmov_divergent": 150, "ecoli": 60, "unplaceable": 3}

READ_LEN = 150
INSERT_MEAN, INSERT_SD = 350, 50

#: Measured on the discriminating dataset, reference-free from self-overlapping pairs
#: (results/berlin_quality/, CLAUDE.md). Not an assumption and not a Q-string estimate,
#: which is an upper bound conservative by about 2.6x.
ERR_BODY, ERR_LAST = 0.00068, 0.0019

#: The nucleotide substitution rate applied to PMMoV to make the divergent arm.
#:
#: Chosen so that the resulting PROTEIN identity to the real PMMoV lands inside the band
#: where a protein search still finds a relative and a nucleotide search mostly does not
#: - the recognition wall this project has measured is at protein level around 10-20%
#: divergence (SPEC 4.2), and the holdout positive control sits at 66.5% protein identity
#: over its replicase. The rate is the input; the identity is MEASURED after mutating and
#: written into the manifest, because the mapping from one to the other runs through the
#: genetic code and the rejection rule below and is not worth predicting.
#: 0.15 was chosen by measuring: 0.22 put the arm at 56% protein identity, below the
#: holdout's 66.5%, and the holdout is the one distance in this repo that nobody picked.
DIVERGENT_NT_RATE = 0.15

STOPS = {"TAA", "TAG", "TGA"}
EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=nuccore&id={acc}&rettype={rettype}&retmode=text")

CODON_TABLE = {}
for _i, _b1 in enumerate("TCAG"):
    for _b2 in "TCAG":
        for _b3 in "TCAG":
            CODON_TABLE[_b1 + _b2 + _b3] = (
                "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
            )[("TCAG".index(_b1) * 16 + "TCAG".index(_b2) * 4 + "TCAG".index(_b3))]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def translate(seq: str) -> str:
    return "".join(CODON_TABLE.get(seq[i:i + 3], "X") for i in range(0, len(seq) - 2, 3))


def read_fasta(text: str) -> list[tuple[str, str]]:
    out, name, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:], []
        elif name is not None:
            buf.append(line.strip().upper())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 70) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")


def fetch(acc: str, rettype: str, refs: Path | None) -> str:
    """Fetch one record, from a local directory if one was given.

    `--refs` exists so that a rebuild does not depend on NCBI being reachable, and so
    that the sandbox this repo is usually worked in - which cannot reach NCBI at all
    (CLAUDE.md) - can still rebuild the fixture from downloaded copies.
    """
    if refs is not None:
        p = refs / (f"{acc}_aa.fa" if rettype == "fasta_cds_aa"
                    else f"{acc}_cds.fa" if rettype == "fasta_cds_na" else f"{acc}.fa")
        if p.exists():
            return p.read_text()
    url = EFETCH.format(acc=acc, rettype=rettype)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                text = r.read().decode()
            if text.startswith(">"):
                time.sleep(0.4)
                return text
        except Exception:
            pass
        time.sleep(2 ** attempt)
    raise SystemExit(f"could not fetch {acc} {rettype}; pass --refs with local copies")


def cds_spans(genome: str, cds_na: str) -> list[tuple[int, int, str]]:
    """Locate each annotated CDS in the genome by searching for its own sequence.

    The CDS records are fetched as nucleotide sequence rather than as coordinates because
    the FASTA headers NCBI returns for `fasta_cds_na` carry no `[location=...]` on these
    records. A CDS that is spliced or wraps the origin will not be found as one substring,
    and this RAISES in that case rather than silently leaving the ORF unprotected: the
    rejection rule below is what keeps the divergent arm free of introduced stop codons,
    and an ORF it does not know about would not be protected.
    """
    spans = []
    for name, seq in read_fasta(cds_na):
        i = genome.find(seq)
        if i < 0:
            raise SystemExit(f"CDS {name.split()[0]} is not a contiguous substring of its "
                             f"genome; this fixture's mutation step cannot protect it")
        spans.append((i, i + len(seq), name.split()[0]))
    return spans


def diverge(genome: str, spans: list[tuple[int, int, str]], rate: float,
            rng: random.Random) -> str:
    """Substitute nucleotides at `rate`, never creating a stop codon or losing a start.

    WHY THE REJECTION RULE. An unconstrained substitution at 22% puts a premature stop in
    every ORF of a 6 kb genome several times over, and the divergent arm would then
    assemble into a contig whose ORFs are fragments. The point of the arm is a contig that
    a protein search still places, which needs the reading frames intact. So a proposed
    substitution is rejected when it would make a stop codon inside any CDS, or touch the
    first codon of one. Synonymous and non-synonymous changes are otherwise accepted
    alike: biasing towards one would be inventing a selection regime.
    """
    g = list(genome)
    frames = {}
    for start, end, _name in spans:
        for pos in range(start, end):
            frames.setdefault(pos, []).append((start, pos - start))
    order = list(range(len(g)))
    rng.shuffle(order)
    target = int(round(rate * len(g)))
    made = 0
    for pos in order:
        if made >= target:
            break
        for base in rng.sample("ACGT", 4):
            if base == g[pos]:
                continue
            ok = True
            for start, offset in frames.get(pos, []):
                if offset < 3:                      # the start codon stays a start codon
                    ok = False
                    break
                c0 = start + (offset // 3) * 3
                codon = g[c0:c0 + 3].copy()
                codon[offset % 3] = base
                if "".join(codon) in STOPS:
                    ok = False
                    break
            if ok:
                g[pos] = base
                made += 1
                break
    return "".join(g)


def identity(a: str, b: str) -> float:
    n = min(len(a), len(b))
    return sum(1 for i in range(n) if a[i] == b[i]) / n if n else 0.0


def simulate(records: list[tuple[str, str, int]], rng: random.Random
             ) -> tuple[list[tuple[str, str, str]], list[tuple[str, str, str]]]:
    """Simulate paired reads at the stated per-genome coverage.

    Reads are drawn from both strands, which matters because a library drawn from one
    strand lets an assembler off a job it has to do on real data, and the repo has an
    audited strand-aware labeller for exactly that reason (results/labeller_audit/).

    The read NAME CARRIES THE SOURCE GENOME. That is only usable by a test - it is ground
    truth, which never reaches a report (SPEC 13.4) - and it is what lets the end-to-end
    test say which class a missing section lost, rather than only that the section is
    missing.
    """
    r1, r2 = [], []
    comp = str.maketrans("ACGTN", "TGCAN")
    qual_body = chr(33 + 37)
    for label, seq, cov in records:
        n_pairs = max(1, int(round(cov * len(seq) / (2 * READ_LEN))))
        for i in range(n_pairs):
            for _attempt in range(20):
                ins = int(round(rng.gauss(INSERT_MEAN, INSERT_SD)))
                if 2 * READ_LEN <= ins <= len(seq):
                    break
            else:
                ins = min(len(seq), 2 * READ_LEN)
            start = rng.randrange(0, len(seq) - ins + 1)
            frag = seq[start:start + ins]
            if rng.random() < 0.5:
                frag = frag.translate(comp)[::-1]
            a = frag[:READ_LEN]
            b = frag[-READ_LEN:].translate(comp)[::-1]
            name = f"{label}_{i}"
            r1.append((name, errors(a, rng), qual_body * len(a)))
            r2.append((name, errors(b, rng), qual_body * len(b)))
    order = list(range(len(r1)))
    rng.shuffle(order)
    return [r1[i] for i in order], [r2[i] for i in order]


def errors(seq: str, rng: random.Random) -> str:
    out = list(seq)
    for i, base in enumerate(out):
        p = ERR_LAST if i == len(out) - 1 else ERR_BODY
        if rng.random() < p:
            out[i] = rng.choice([b for b in "ACGT" if b != base])
    return "".join(out)


def write_fastq(path: Path, reads: list[tuple[str, str, str]], mate: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0 so that a rebuild is byte-identical: gzip stamps the wall clock into its
    # header otherwise, and a fixture whose checksum changes on every rebuild cannot be
    # pinned by the manifest that exists to pin it.
    # THE MATE GOES AFTER THE SPACE, NOT INTO THE NAME. Both mates of a pair carry the
    # same identifier in real Illumina output, and every tool downstream takes the token
    # before the first space as the read's name. Written as `name/1`, the mate suffix
    # became part of the name - and Kraken2, which strips a trailing /1 or /2 when it is
    # given paired input, then returned a set of names that matched nothing else in the
    # run. The cascade's own accounting caught it, with a removed set larger than the
    # block it was removing from and a negative residue, which is the impossible case an
    # assertion is for.
    payload = "".join(f"@{name} {mate}:N:0:1\n{seq}\n+\n{qual}\n"
                      for name, seq, qual in reads)
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=9,
                                                mtime=0) as fh:
        fh.write(payload.encode())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--refs", type=Path, default=None,
                    help="directory of already-downloaded records, used in preference to "
                         "NCBI so a rebuild works offline")
    ap.add_argument("--out", type=Path, default=HERE)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    out = args.out

    genomes, proteins, cds = {}, {}, {}
    for key, meta in ACCESSIONS.items():
        acc = meta["accession"]
        genomes[key] = read_fasta(fetch(acc, "fasta", args.refs))[0][1]
        if key != "ecoli":
            proteins[key] = read_fasta(fetch(acc, "fasta_cds_aa", args.refs))
            cds[key] = fetch(acc, "fasta_cds_na", args.refs)

    lo, hi = ECOLI_SLICE
    ecoli = genomes["ecoli"][lo - 1:hi]

    spans = cds_spans(genomes["pmmov"], cds["pmmov"])
    divergent = diverge(genomes["pmmov"], spans, DIVERGENT_NT_RATE, rng)

    # MEASURED, not stated: per-ORF protein identity of the divergent arm against the
    # real PMMoV proteins that the toy protein index carries. These are the numbers the
    # end-to-end test reads to know what the divergent section should be able to show,
    # and they are a property of the mutation above, so they are recomputed on every
    # rebuild rather than copied forward.
    measured = []
    for start, end, name in spans:
        real = translate(genomes["pmmov"][start:end]).rstrip("*")
        mut = translate(divergent[start:end]).rstrip("*")
        # STOPS ARE COMPARED, NOT COUNTED. PMMoV's replicase is a readthrough product, so
        # the shorter 126K CDS ends at a stop codon that sits in frame INSIDE the longer
        # ORF's span: the real ORF 1 translation already carries an internal `*`. An
        # assertion that the mutated ORF has no internal stop therefore fails on correct
        # output, which is what it did on the first build. What has to hold is that the
        # mutation introduced none, so the stop POSITIONS are compared against the real
        # translation's own.
        real_stops = {i for i, c in enumerate(real) if c == "*"}
        mut_stops = {i for i, c in enumerate(mut) if c == "*"}
        measured.append({"orf": name, "aa": len(real),
                         "protein_identity": round(identity(real, mut), 4),
                         "nucleotide_identity": round(
                             identity(genomes["pmmov"][start:end], divergent[start:end]), 4),
                         "internal_stops_in_real": sorted(real_stops),
                         "stops_introduced": sorted(mut_stops - real_stops)})
    assert not any(m["stops_introduced"] for m in measured), measured

    # Sequence in no index at all. Drawn from the fixture's own nucleotide composition
    # rather than uniformly, so that a composition field cannot separate it from the real
    # arms for free - composition as a novelty signal is dead in this repo (SPEC 11) and
    # a fixture that revives it by accident would be misleading.
    pool = genomes["ms2"] + ecoli
    freq = [pool.count(b) for b in "ACGT"]
    unplaceable = "".join(rng.choices("ACGT", weights=freq, k=4000))

    write_fasta(out / "genomes/known_ms2.fa",
                [(f"{ACCESSIONS['ms2']['accession']} MS2, known", genomes["ms2"])])
    write_fasta(out / "genomes/divergent_pmmov.fa",
                [("fixture_divergent_pmmov "
                  f"mutated from {ACCESSIONS['pmmov']['accession']} at "
                  f"nt_rate={DIVERGENT_NT_RATE}", divergent)])
    write_fasta(out / "genomes/background_ecoli.fa",
                [(f"{ACCESSIONS['ecoli']['accession']}:{lo}-{hi} E. coli rrnB window",
                  ecoli)])
    write_fasta(out / "genomes/unplaceable.fa",
                [("fixture_unplaceable composition-matched, in no index", unplaceable)])

    # The toy indexes, as source sequence. They are BUILT in CI by build_dbs.sh, not
    # committed as binaries: a .dmnd or a Kraken2 hash is a build artefact of one tool
    # version and the pins live in TOOLS.md.
    write_fasta(out / "db/nt_known.fa",
                [(f"{ACCESSIONS['ms2']['accession']} MS2", genomes["ms2"]),
                 (f"{ACCESSIONS['ecoli']['accession']}:{lo}-{hi} E. coli", ecoli)])
    # The protein index carries the REAL PMMoV, never the mutated one: that is what makes
    # the divergent arm divergent rather than present.
    write_fasta(out / "db/prot_known.faa",
                [(n.split()[0], s) for n, s in proteins["ms2"] + proteins["pmmov"]])
    # The toy SILVA is the rRNA operon window only, which is what SILVA is: an rRNA
    # index, not a genome index.
    write_fasta(out / "db/silva_toy.fa",
                [(f"{ACCESSIONS['ecoli']['accession']}:{lo}-{hi} rrnB window", ecoli)])

    arms = [("ms2", genomes["ms2"], COVERAGE["ms2"]),
            ("pmmovdiv", divergent, COVERAGE["pmmov_divergent"]),
            ("ecoli", ecoli, COVERAGE["ecoli"]),
            ("unplaceable", unplaceable, COVERAGE["unplaceable"])]
    r1, r2 = simulate(arms, rng)
    write_fastq(out / "reads/fixture_R1.fastq.gz", r1, 1)
    write_fastq(out / "reads/fixture_R2.fastq.gz", r2, 2)

    files = sorted(p for p in out.rglob("*")
                   if p.is_file() and p.name not in ("manifest.json", "make_fixture.py")
                   and "__pycache__" not in p.parts)
    manifest = {
        "what_this_is": "the end-to-end fixture: a ten-thousand-pair library carrying all "
                        "three of the report's output classes, with toy indexes built from "
                        "its own sequence. It is a test of whether every report section is "
                        "filled, and it is NOT a benchmark: no number from it belongs in "
                        "results/.",
        "script": "tests/fixtures/make_fixture.py",
        "script_sha256": sha256(Path(__file__)),
        "seed": args.seed,
        "built_on": date.today().isoformat(),
        "accessions": ACCESSIONS,
        "ecoli_slice_1based_inclusive": list(ECOLI_SLICE),
        "coverage": COVERAGE,
        "read_length": READ_LEN,
        "insert_mean_sd": [INSERT_MEAN, INSERT_SD],
        "error_rate": {"body": ERR_BODY, "last_base": ERR_LAST,
                       "basis": "measured reference-free from self-overlapping pairs of "
                                "SRR23973193 (results/berlin_quality/, CLAUDE.md)"},
        "divergent_arm": {
            "nucleotide_substitution_rate_applied": DIVERGENT_NT_RATE,
            "rejection_rule": "a substitution that would make a stop codon inside any "
                              "annotated CDS, or touch a start codon, is rejected",
            "measured_per_orf": measured,
            "basis": "measured on this build, not stated: the mapping from a nucleotide "
                     "rate to a protein identity runs through the genetic code and the "
                     "rejection rule.",
        },
        "pairs": len(r1),
        "files": {str(p.relative_to(out)): {"bytes": p.stat().st_size,
                                            "sha256_16": sha256(p)[:16]} for p in files},
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"{len(r1)} pairs; divergent arm protein identity per ORF: "
          f"{[m['protein_identity'] for m in measured]}")
    for p in files:
        print(f"  {p.relative_to(out)}  {p.stat().st_size:,} bytes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
