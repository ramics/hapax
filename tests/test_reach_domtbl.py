"""HMMER's two coordinate pairs, and which of them is our ORF.

WHY THIS TEST EXISTS, with the bug it is a regression for. `--domtblout` reports an
alignment's "hmm coord" pair on the QUERY and its "ali coord" pair on the TARGET. Which of
those is our ORF depends on the program: `hmmsearch` searches a profile library against our
sequences, so the ORF is the target; `jackhmmer` searches our sequence against a sequence
database, so the ORF is the query. One mapping was used for both on run
`i-0672a2ca4a782aa21`, which put the Pfam profile's length into the ORF's `qlen` column and
the subject's alignment span into the ORF's. The only visible symptom was a query coverage
of 1.156 - a number above 1.0 in a table of thousands of rows.

So the fixtures here are real lines from that run's own `pfam_ga.domtbl`, and the assertions
are the two things that would have caught it: coverage can never exceed 1.0, and `qlen` has
to be the ORF's length rather than the profile's. `leak__orf4` is the useful row - a 1,132
residue ORF matching a 173-position profile over 200 residues - because every quantity in it
is different from every other, so a swap cannot pass by coincidence. A row where the ORF and
the profile happen to be similar lengths would not have failed.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]

#: Verbatim from `results/reach/` run i-0672a2ca4a782aa21, `pfam_ga.domtbl`. hmmsearch
#: against Pfam-A: target is the ORF, query is the profile.
HMMSEARCH_LINES = """\
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
leak__orf4           -           1132 Chropara_Vmeth       PF19223.7    173   5.4e-26   83.7   2.6   1   1   7.9e-28   7.1e-26   83.3   0.4     3   172    94   293    92   294 0.90 -
ctl6__orf3           -            150 Levi_coat            PF01819.25   135   2.9e-46  149.4   1.9   1   1   3.5e-47   3.1e-46  149.3   1.9     1   130     2   145     2   150 0.93 -
"""

#: The same format the other way round. jackhmmer against the sequence set: query is the
#: ORF, target is a database protein.
JACKHMMER_LINES = """\
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
YP_010769254.1       -            552 div1__orf1           -            555  4.0e-171  576.6   0.0   1   1  1.0e-172  4.5e-171  575.0   0.0     6   550     4   548     2   551 0.95 -
"""


@pytest.fixture(scope="module")
def arms():
    spec = importlib.util.spec_from_file_location(
        "reach_arms", ROOT / "experiments" / "reach_arms.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["reach_arms"] = module
    spec.loader.exec_module(module)
    return module


def test_hmmsearch_reports_the_orfs_own_length_and_span(arms, tmp_path):
    """The ORF is the target, so its length is `tlen` and its span is the ali coords."""
    path = tmp_path / "pfam_ga.domtbl"
    path.write_text(HMMSEARCH_LINES)
    rows = arms.parse_domtbl(path, arm="hmmsearch_ga", program="hmmsearch --cut_ga",
                             db="Pfam-A", orf_is_target=True, accession_as_subject=True)
    assert len(rows) == 2
    big = next(r for r in rows if r["orf"] == "leak__orf4")

    # The row that cannot pass a swap by coincidence: every quantity differs.
    assert big["qlen"] == 1132, "qlen is the profile's length, not the ORF's"
    assert big["slen"] == 173
    assert (big["qstart"], big["qend"]) == (94, 293), "the ORF's span is the ali coords"
    assert (big["sstart"], big["send"]) == (3, 172), "the profile's span is the hmm coords"
    assert big["aln"] == 200
    assert big["qcov"] == pytest.approx(200 / 1132)
    assert big["subject"] == "PF19223.7", "the accession is the stable identifier"
    assert big["profile_name"] == "Chropara_Vmeth"
    assert big["pident"] is None, "a profile has no pairwise identity to report"


def test_jackhmmer_reports_the_orfs_own_length_and_span(arms, tmp_path):
    """The ORF is the query, so the same two pairs mean the opposite things."""
    path = tmp_path / "jackhmmer.tbl"
    path.write_text(JACKHMMER_LINES)
    rows = arms.parse_domtbl(path, arm="jackhmmer", program="jackhmmer -N 5",
                             db="swissprot+viral", orf_is_target=False)
    assert len(rows) == 1
    row = rows[0]
    assert row["orf"] == "div1__orf1"
    assert row["subject"] == "YP_010769254.1"
    assert row["qlen"] == 555, "qlen is the subject's length, not the ORF's"
    assert row["slen"] == 552
    assert (row["qstart"], row["qend"]) == (6, 550), "the ORF's span is the hmm coords here"
    assert (row["sstart"], row["send"]) == (4, 548)
    assert row["qcov"] == pytest.approx(545 / 555)


@pytest.mark.parametrize("orf_is_target", [True, False])
def test_coverage_never_exceeds_one(arms, tmp_path, orf_is_target):
    """The assertion that would have caught the original bug, both ways round.

    A coverage above 1.0 means the span and the length came from different sequences. It is
    the only cheap invariant this format offers, so it is tested rather than trusted.
    """
    path = tmp_path / f"dom_{orf_is_target}.domtbl"
    path.write_text(HMMSEARCH_LINES if orf_is_target else JACKHMMER_LINES)
    rows = arms.parse_domtbl(path, arm="x", program="x", db="x",
                             orf_is_target=orf_is_target)
    assert rows
    for row in rows:
        assert 0 < row["qcov"] <= 1.0, f"{row['orf']} has coverage {row['qcov']}"
