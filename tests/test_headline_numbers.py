"""Every headline number, recomputed from the committed files it came from.

WHAT THIS IS, AND THE WEAKER CLAIM IT MAKES. PLAN.md's last definition-of-done item is a
regression test that re-running the Berlin and Yeongdok inputs reproduces their committed
headline numbers. Re-running them is tens of hours on 8.7 M and 32.1 M pairs and needs
indexes this container cannot hold, so what runs here is the half that can: each headline
is RECOMPUTED FROM THE COMMITTED PER-SEED FILES, and the prose that quotes it is checked
to still say the same thing. That does not prove the pipeline reproduces. It does catch
the failure this repo has actually had, twice: prose drifting from the outputs under it,
which is why figures are generated from committed files and never hand-written. The
stronger claim needs a run, and `results/PROVENANCE.md` is where it is tracked.

WHY EACH NUMBER NAMES ITS OWN FILE. The same quantity has two values here. The cascade's
union reads 77.10% over `results/subtract_recut/` and 89.18% over
`results/subtract_known/`, because the earlier run subtracted contigs it could not
identify and swallowed tier 2; the recut is the headline and the inflated figure is kept
because what it measures is the mistake. A test that averaged "the cascade's seeds" would
mix them. So a claim here is a (number, file, prose location) triple, and the file is part
of the claim.

THE SPREAD IS AN OBSERVED RANGE, NOT A DISPERSION ESTIMATE. The min and max asserted
below are the five seeds' own, and five seeds give an honest range and a useless variance:
the sampling error of a five-sample standard deviation swamps what it would measure, and
off five of this library's blocks the same correction reverses sign. So this asserts what
was seen and derives nothing from it.
"""
from __future__ import annotations

import json
import statistics
from dataclasses import dataclass
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _cascade(paths: list[Path]) -> dict:
    """Mean and observed range of the cascade's own per-seed fields.

    Percentages are of READ PAIRS, the occurrence-weighted unit, which is what each file
    records: averaging over distinct k-mers instead differs by two orders of magnitude
    for conserved sequence, and a read sees the occurrence-weighted one.
    """
    seeds = [json.loads(p.read_text()) for p in paths]
    assert seeds, "no committed seed files"
    union = [s["union_removed"] for s in seeds]
    residue = [s["residue_pairs"] / s["pairs"] for s in seeds]
    out = {"seeds": len(seeds),
           "union": 100 * statistics.mean(union),
           "union_low": 100 * min(union), "union_high": 100 * max(union),
           "residue": 100 * statistics.mean(residue)}
    if all("tier3_pairs" in s for s in seeds):
        out["tier3"] = 100 * statistics.mean(s["tier3_pairs"] / s["pairs"]
                                             for s in seeds)
    return out


@dataclass(frozen=True)
class Claim:
    """One headline: what it says, where the number comes from, where it is quoted."""
    label: str
    #: The committed files it is recomputed from, relative to the repo root.
    files: tuple[str, ...]
    #: field -> the value the prose states, to two decimals.
    expect: dict
    #: The exact strings that must still appear in the prose, by file. A number quoted
    #: with no entry here can drift without this test noticing, which is why the prose
    #: side is asserted rather than assumed.
    quoted: dict


CLAIMS = (
    Claim(
        label="Berlin SRR23973193: identification removes 77.10%, residue 22.90%",
        files=tuple(f"results/subtract_recut/subtract_s{s}.json" for s in range(1, 6)),
        expect={"seeds": 5, "union": 77.10, "union_low": 77.02, "union_high": 77.28,
                "residue": 22.90},
        quoted={"SPEC.md": ["**Union 77.10% (77.02–77.28 across seeds). Residue 22.90%.**"],
                "results/README.md": ["**77.10%**"]},
    ),
    Claim(
        label="the pre-recut run that subtracted unidentifiable contigs: 89.18%",
        files=tuple(f"results/subtract_known/subtract_s{s}.json" for s in range(1, 6)),
        expect={"seeds": 5, "union": 89.18},
        quoted={"results/subtract_known/REPORT.md": ["89.18%"]},
    ),
    Claim(
        label="Yeongdok SRR23702504: 78.42%, residue 21.58%, tier 3 11.74%",
        files=tuple(f"results/prod_SRR23702504/run_s{s}/tiers.json" for s in range(1, 6)),
        expect={"seeds": 5, "union": 78.42, "union_low": 78.31, "union_high": 78.59,
                "residue": 21.58, "tier3": 11.74},
        quoted={"results/README.md": ["**78.42%**", "21.58%", "11.74%"]},
    ),
)


@pytest.mark.parametrize("claim", CLAIMS, ids=[c.label for c in CLAIMS])
def test_the_headline_is_what_its_committed_files_say(claim):
    paths = [ROOT / f for f in claim.files]
    missing = [p for p in paths if not p.exists()]
    if missing:
        pytest.skip(f"not in this tree: {[str(p) for p in missing]}")
    got = _cascade(paths)
    for field, want in claim.expect.items():
        assert field in got, f"{claim.label}: {field} is not computed from these files"
        if field == "seeds":
            assert got[field] == want, (
                f"{claim.label}: {got[field]} seeds committed, not {want}. A mixture is "
                f"not a five-seed result")
            continue
        assert round(got[field], 2) == want, (
            f"{claim.label}: {field} recomputes to {got[field]:.4f} and the prose says "
            f"{want}. Either the files changed or the prose drifted from them")


@pytest.mark.parametrize("claim", CLAIMS, ids=[c.label for c in CLAIMS])
def test_the_prose_still_quotes_the_number(claim):
    """The other direction. A result file and a sentence about it can disagree either way,
    and this repo has had to withdraw claims whose prose drifted from its outputs."""
    for relative, strings in claim.quoted.items():
        path = ROOT / relative
        if not path.exists():
            pytest.skip(f"{relative} is not in this tree")
        text = path.read_text()
        for needle in strings:
            assert needle in text, (
                f"{relative} no longer says {needle!r}, which {claim.label} rests on. If "
                f"the number moved, move it in the file the test reads too")


def test_the_discriminating_library_facts_are_the_ones_every_number_was_measured_on():
    """The library's own facts, which are not derivable from any result file.

    A trim and a block size decide which bases and which reads a number was measured
    over, so they belong with the numbers. They live in `experiments/libraries.json` with
    how each was established; this asserts the values rather than assuming them, which is
    the same reason `tests/test_library_facts.py` exists.
    """
    libs = json.loads((ROOT / "experiments" / "libraries.json").read_text())["libraries"]
    berlin = libs["SRR23973193"]
    assert berlin["pairs"] == 8_731_224
    assert berlin["trim_r2"] == 14
    assert berlin["block"] == 1_500_000
    # A second operating point of the same quantity gets its own key rather than a
    # reconciliation, and each says so: 16 for Kraken2 and 16 for mapback, the latter
    # because that is what `results/divergent_known/` was produced with.
    assert berlin["trim_r2_for_kraken2"] == 16
    assert berlin["trim_r2_for_mapback"] == 16
    for key in ("trim_r2", "block", "trim_r2_for_kraken2", "trim_r2_for_mapback"):
        assert berlin[f"{key}_basis"], (
            f"{key} is stated for SRR23973193 with no account of where it came from, "
            f"which is the default this repo removed from twenty scripts")
