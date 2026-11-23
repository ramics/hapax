"""The codon arm's walk, against alignments whose answer is known by construction.

WHY THIS TEST EXISTS. `arm_codon` reconstructs a codon alignment from three things that do
not agree about their units: a protein alignment with gaps in it, the query's CDS in
nucleotides, and the subject's coordinates in a genome - which can run backwards. An
off-by-one in any of them produces dN/dS values that are wrong and completely plausible,
on a box with no SSH, in a table nobody can check by eye. So the walk is tested here,
where the synonymous and non-synonymous differences were put in on purpose and can be
counted by hand.

The four cases are the four things that can go wrong:

1. **Identical sequences** - every count must be zero and dS must be undefined rather than
   zero-divided, because Jukes-Cantor of zero divergence is zero and the ratio is then 0/0.
2. **Synonymous differences only** - the third-position changes must land in `syn_diff` and
   nowhere else, which is the test that the reading frame is right. A frame shifted by one
   turns every synonymous change into a non-synonymous one, so this case fails loudly on
   exactly the error that is otherwise invisible.
3. **A subject on the reverse strand** - the walk steps backwards and reverse-complements,
   and getting that wrong gives a garbage codon that the stop filter mostly swallows. Here
   the reverse-strand case is the same alignment as case 2 with the genome flipped, so it
   must give the identical answer.
4. **Saturation** - a heavily diverged pair where the Jukes-Cantor argument goes
   non-positive, which must come back as `None` and `synonymous_saturated=True` rather than
   as a finite number. That behaviour is the whole point of the arm: a clipped dS would
   turn "this cannot be measured here" into a measurement.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

TAB_FIELDS = ("qseqid", "sseqid", "pident", "length", "qlen", "slen", "qstart", "qend",
              "sstart", "send", "evalue", "bitscore", "qcovhsp", "qseq", "sseq", "sframe")

CODE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L",
    "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V",
    "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P",
    "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "CAT": "H",
    "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R",
    "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def translate(nt: str) -> str:
    return "".join(CODE[nt[i:i + 3]] for i in range(0, len(nt) - len(nt) % 3, 3))


def revcomp(nt: str) -> str:
    return nt[::-1].translate(str.maketrans("ACGT", "TGCA"))


@pytest.fixture(scope="module")
def arms():
    spec = importlib.util.spec_from_file_location(
        "reach_arms", ROOT / "experiments" / "reach_arms.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["reach_arms"] = module
    spec.loader.exec_module(module)
    return module


def run_codon(arms, tmp_path: Path, query_cds: str, genome: str, *,
              sstart: int, send: int, translation_table: int = 11) -> dict:
    """One alignment through `arm_codon`, with the tblastn table written by hand.

    The protein strings are translated from the two nucleotide sequences rather than typed,
    so the test cannot disagree with its own fixtures - which would make a failure here
    ambiguous between the walk and the test.
    """
    q_aa = translate(query_cds)
    subject_cds = genome[sstart - 1:send] if send >= sstart else revcomp(
        genome[send - 1:sstart])
    s_aa = translate(subject_cds)
    assert len(q_aa) == len(s_aa), "the fixture's two sequences are not the same length"

    out = tmp_path
    out.mkdir(parents=True, exist_ok=True)
    (out / "orfs.json").write_text(json.dumps([dict(
        orf="q1", query="c1", arm="divergent", pile="real", start=1,
        end=len(query_cds), strand=1, aa_len=len(q_aa),
        translation_table=translation_table, partial=False, cds_nt=query_cds)]))
    row = dict(qseqid="q1", sseqid="g1", pident=100.0 * sum(
        1 for a, b in zip(q_aa, s_aa) if a == b) / len(q_aa), length=len(q_aa),
        qlen=len(q_aa), slen=len(genome), qstart=1, qend=len(q_aa), sstart=sstart,
        send=send, evalue=1e-40, bitscore=200.0, qcovhsp=100.0, qseq=q_aa, sseq=s_aa,
        sframe=1 if send >= sstart else -1)
    (out / "tblastn.tsv").write_text(
        "\t".join(str(row[f]) for f in TAB_FIELDS) + "\n")

    rows = arms.arm_codon({}, dict(viral_nt=out / "viral_nt",
                                   viral_seqs={"g1": genome}), out, 1)
    assert len(rows) == 1, f"expected one row, got {len(rows)}"
    return rows[0]


def test_identical_sequences_have_no_differences_and_an_undefined_ratio(arms, tmp_path):
    """Zero divergence: no differences anywhere, and dS of zero rather than a ratio."""
    cds = "ATGGCTAAAGGTTTCCCAGACGAATGGCGTCATACTGTTAAC"
    genome = "GGGCCC" + cds + "TTTAAA"
    row = run_codon(arms, tmp_path, cds, genome, sstart=7, send=6 + len(cds))
    assert row["codon_pairs"] == len(cds) // 3
    assert row["p_syn"] == 0.0 and row["p_nonsyn"] == 0.0
    assert row["d_syn"] == 0.0 and row["d_nonsyn"] == 0.0
    # Jukes-Cantor of zero is zero, which is defined; the ratio is what is not, and the
    # field must be None rather than a ZeroDivisionError or an invented number.
    assert row["dnds"] is None
    assert row["synonymous_saturated"] is False


def test_third_position_changes_are_all_synonymous(arms, tmp_path):
    """The frame test. A frame off by one makes every one of these non-synonymous."""
    cds = "ATGGCTAAAGGTTTCCCAGACGAATGGCGTCATACTGTTAAC"
    # Four fourfold-degenerate codons changed at the third position only: GCT->GCC (Ala),
    # GGT->GGC (Gly), CGT->CGC (Arg), ACT->ACC (Thr). Chosen because all four stay the same
    # residue under every table this arm uses.
    subject = list(cds)
    for pos, base in ((5, "C"), (11, "C"), (29, "C"), (35, "C")):
        subject[pos] = base
    subject = "".join(subject)
    assert translate(subject) == translate(cds), "the fixture changed a residue"
    genome = "AC" + subject + "GT"
    row = run_codon(arms, tmp_path, cds, genome, sstart=3, send=2 + len(subject))
    assert row["p_nonsyn"] == 0.0, "a synonymous-only pair produced non-synonymous changes"
    assert row["p_syn"] > 0.0
    assert row["d_nonsyn"] == 0.0
    assert row["dnds"] == 0.0


def test_reverse_strand_subject_gives_the_same_answer(arms, tmp_path):
    """Same alignment, genome flipped. The walk must step backwards and complement."""
    cds = "ATGGCTAAAGGTTTCCCAGACGAATGGCGTCATACTGTTAAC"
    subject = list(cds)
    for pos, base in ((5, "C"), (11, "C"), (29, "C"), (35, "C")):
        subject[pos] = base
    subject = "".join(subject)

    forward_genome = "AC" + subject + "GT"
    forward = run_codon(arms, tmp_path / "fwd", cds, forward_genome,
                        sstart=3, send=2 + len(subject))

    # The same subject CDS on the other strand: the genome is the reverse complement, and
    # the alignment's subject coordinates therefore run downwards.
    reverse_genome = revcomp(forward_genome)
    start = len(reverse_genome) - 2            # where the CDS begins, counting from 1
    row = run_codon(arms, tmp_path / "rev", cds, reverse_genome,
                    sstart=start, send=start - len(subject) + 1)

    assert row["codon_pairs"] == forward["codon_pairs"]
    assert row["p_syn"] == forward["p_syn"]
    assert row["p_nonsyn"] == forward["p_nonsyn"] == 0.0


def test_saturated_synonymous_sites_come_back_undefined(arms, tmp_path):
    """The behaviour the arm exists for: unmeasurable is reported, not clipped.

    Jukes-Cantor's argument is 1 - 4p/3, so a synonymous proportion at or above 0.75 has no
    correction. A clipped dS there would turn saturation into a finite number and hide
    exactly the thing this arm was built to test.
    """
    # Every codon fourfold degenerate, every third position changed to a different base:
    # the synonymous proportion goes to its ceiling while no residue changes at all.
    codons = ["GCT", "GGT", "CGT", "ACT", "CCT", "GTT", "TCT", "CTT"] * 4
    cds = "ATG" + "".join(codons)
    subject_codons = [c[:2] + ("C" if c[2] != "C" else "A") for c in codons]
    subject = "ATG" + "".join(subject_codons)
    assert translate(subject) == translate(cds)
    row = run_codon(arms, tmp_path, cds, subject, sstart=1, send=len(subject))
    assert row["p_nonsyn"] == 0.0
    assert row["p_syn"] >= 0.75, f"fixture is not saturated: p_syn={row['p_syn']}"
    assert row["d_syn"] is None, "a saturated dS was clipped to a finite value"
    assert row["synonymous_saturated"] is True
    assert row["dnds"] is None
