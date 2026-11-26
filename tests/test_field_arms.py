"""Every field's three arms, and the two that were computed and never aggregated.

A field ships with a null, which says what it reads on sequence carrying no information,
and with a known-positive arm, which says whether it reads differently on sequence the
cascade did identify. A field that reads the same on both carries no information about the
residue, and that is its calibration rather than a reason to hide it.

Three fields on a published page read `cannot discriminate here`, and none of the three was
a measurement that failed to separate:

- `gc` had a null identical to itself to every digit, because a composition-matched shuffle
  permutes bases and cannot change GC. That is no null rather than an uninformative one,
  and the two must not render the same.
- `closest_viral` and `closest_sprot` had no null at all. The search ran over the residue
  and over the subtracted reads and never over the shuffled copies, which the same loop was
  already building for the composition fields.
- and every composition field's subtracted-read arm was computed per read and dropped at
  the summary, although the per-read report's own header has promised it since it was
  written.

So: the null arm comes from the shared shuffle and validates its own output; the composition
arms are computed by one accessor per field for all three piles; and the driver's state
column tells "no null exists for this field" apart from "cannot discriminate".
"""
from __future__ import annotations

import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))

import residue_report as rr  # noqa: E402
import shuffle_control  # noqa: E402


def sequence(seed: int, n: int = 150) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def row(gc: float, entropy: float, orf: int, frame_ok: int, dup: float) -> dict:
    """A per-read row of the shape `_row` writes, with only the fields the arms read."""
    return {"gc": gc, "entropy3": entropy, "entropy3_shuffled": entropy - 0.01,
            "frame": {"longest_orf_aa": orf, "frame_identifiable": frame_ok},
            "frame_shuffled": {"longest_orf_aa": orf - 5, "frame_identifiable": 0},
            "repetitive": {"kmer16_duplicate_frac": dup},
            "repetitive_shuffled": {"kmer16_duplicate_frac": 0.0}}


# ---------------------------------------------------------------- the search null arm

def test_the_search_null_is_a_validated_dinucleotide_shuffle_of_the_same_reads():
    pairs = {f"p{i}": (sequence(i), sequence(100 + i)) for i in range(8)}
    null, census = rr.search_null(pairs, random.Random(1))

    assert census["kind"] == "dinuc", "a mononucleotide null is the easier one to beat"
    assert census["offered_pairs"] == 8 and census["shuffled_pairs"] == len(null)
    for name, (s1, s2) in null.items():
        real1, real2 = pairs[name]
        # The invariant a null is published on, checked on the sequence rather than
        # trusted from the code: same length, same dinucleotide counts, not the input.
        shuffle_control.check_shuffle(real1, s1, "dinuc")
        shuffle_control.check_shuffle(real2, s2, "dinuc")


def test_a_pair_the_shuffle_cannot_take_is_counted_and_never_passed_through():
    """A shuffle that hands back its input is real sequence wearing the label of a null.

    2.85% of contigs falling through a broken walk once supplied 304 of 304 of a null's
    significant protein hits. A homopolymer is genuinely invariant, so it is dropped and
    counted here rather than admitted as a shuffle of itself.
    """
    pairs = {"flat": ("A" * 150, "A" * 150), "real": (sequence(3), sequence(4))}
    null, census = rr.search_null(pairs, random.Random(1))

    assert "flat" not in null and "real" in null
    assert census["dropped_pairs"] == 1
    assert (census["mates_no_eulerian_walk"] +
            census["mates_shuffle_returned_input"]) >= 1


# ------------------------------------------------------- the composition field arms

def test_every_composition_field_carries_the_subtracted_read_arm():
    residue = [row(0.50, 0.90, 80, 1, 0.01), row(0.45, 0.88, 70, 0, 0.0)]
    known = [row(0.52, 0.91, 90, 1, 0.02)]

    arms = rr._composition_arms(residue, known)

    for name in ("gc", "entropy3", "longest_orf_aa", "frame_identifiable_frac",
                 "kmer16_duplicate_frac"):
        assert arms[name]["subtracted_reads"] is not None, name
        assert "subtracted_reads_absent_because" not in arms[name], name


def test_gc_reports_no_null_rather_than_a_null_equal_to_itself():
    arms = rr._composition_arms([row(0.5, 0.9, 80, 1, 0.0)], [row(0.5, 0.9, 80, 1, 0.0)])
    assert arms["gc"]["shuffled_null"] is None
    assert "cannot change GC" in arms["gc"]["shuffled_null_absent_because"]
    # and a field a shuffle CAN move still carries its null
    assert arms["entropy3"]["shuffled_null"] is not None


def test_an_absent_known_pile_is_a_stated_absence_and_not_a_zero():
    arms = rr._composition_arms([row(0.5, 0.9, 80, 1, 0.0)], [])
    for name, v in arms.items():
        assert "subtracted_reads" not in v, name
        assert "no known-positive arm" in v["subtracted_reads_absent_because"], name


# ------------------------------------------------------------ the state the page reads

def test_the_state_column_separates_no_null_from_cannot_discriminate(tmp_path):
    """The driver's own reading of those keys, which is what the page prints."""
    import json

    from hapax.driver import Run

    rundir = tmp_path / "run"
    rundir.mkdir()
    (rundir / "fields.json").write_text(json.dumps({"fields": {
        "gc": {"residue": {"mean": 0.5, "sd": 0.1}, "shuffled_null": None,
               "shuffled_null_absent_because": "a shuffle cannot change GC",
               "subtracted_reads": {"mean": 0.52, "sd": 0.1}},
        "closest_viral": {"any_hit_frac": 0.9, "shuffled_null_any_hit_frac": 0.01,
                          "subtracted_reads_any_hit_frac": 0.95},
        "orphan": {"residue": {"mean": 1.0, "sd": 0.1}},
    }}))
    (rundir / "reads.jsonl").write_text(
        json.dumps({"name": "a", "set": "residue"}) + "\n" +
        json.dumps({"name": "b", "set": "known"}) + "\n" +
        json.dumps({"name": "c", "set": "shuffled"}) + "\n")

    check, known_rows = Run.field_calibration(rundir)

    assert known_rows == 1, "a shuffled row is not a known positive"
    assert "no null exists for this field" in check["gc"]["state"]
    assert "subtracted pile" in check["gc"]["state"]
    assert check["closest_viral"] == "calibrated"
    assert "cannot discriminate" in check["orphan"]
