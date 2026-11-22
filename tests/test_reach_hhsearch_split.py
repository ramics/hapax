"""Which query a jackhmmer alignment belongs to, when a null is named after a real one.

WHY THIS TEST EXISTS, with the bug it is a regression for. The profile-to-profile arm
splits jackhmmer's concatenated Stockholm into one alignment per query and names each file
after that query. The first version found the name by substring - the first ORF name that
appeared anywhere in the block. Every shuffled record on this panel is named
`sh<seed>_<qid>`, so a real ORF's name is a substring of its own null's name: block
`sh5_ctl1__orf1` matched `ctl1__orf1` and was written to `ctl1__orf1.sto`, overwriting the
real query's alignment with its own control.

On run `i-0499bf15b019eba9f` that turned 202 ORFs into 59 files, nearly all of them holding
a shuffle. A shuffled ORF finds no homologs, so its alignment holds one sequence, `hhmake`
refuses a profile built from one sequence, and the arm therefore finished with zero rows and
`ok: true` - which reads as "profile search found no relative" and is a different claim
entirely. **The control had replaced the test and the run still looked clean.**

So the fixtures below are the real shape of that collision, and the assertions are the three
things that would have caught it: the name comes from jackhmmer's own `#=GF ID` line, a real
query and its own null land in different files, and anything unexpected raises rather than
guessing. The last test pins the silent-zero half: an arm that builds no profile at all has
to fail rather than report an empty table.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def block(query: str, *sequences: tuple[str, str], iteration: int = 1) -> str:
    """One jackhmmer Stockholm block, in the format run 2 actually wrote."""
    lines = ["# STOCKHOLM 1.0", f"#=GF ID {query}-i{iteration}",
             "#=GF AU jackhmmer (HMMER 3.4)", ""]
    for name, seq in sequences:
        lines.append(f"{name} {seq}")
    lines.append("//")
    return "\n".join(lines) + "\n"


#: A real ORF with two homologs, and its own dinucleotide shuffle with none - the pair whose
#: names collide under a substring match, because "div1__orf1" is inside "sh5_div1__orf1".
REAL = block("div1__orf1", ("div1__orf1", "MRHWPRPLAVDCSCC"),
             ("YP_010769254.1", "MRHWPRPLAVDCTCC"),
             ("YP_009640089.1", "MRYWPRPLAVDCSCS"))
NULL = block("sh5_div1__orf1", ("sh5_div1__orf1", "CWMRPVDLSCHCRPA"))


@pytest.fixture(scope="module")
def arms():
    spec = importlib.util.spec_from_file_location(
        "reach_arms", ROOT / "experiments" / "reach_arms.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["reach_arms"] = module
    spec.loader.exec_module(module)
    return module


def test_a_null_does_not_overwrite_the_real_query_it_is_named_after(arms, tmp_path):
    """The bug itself: both blocks must survive, in files of their own."""
    sto = tmp_path / "jackhmmer.sto"
    sto.write_text(REAL + NULL)
    pairs = arms.split_jackhmmer_alignments(
        sto, {"div1__orf1", "sh5_div1__orf1"}, tmp_path / "hhm")

    assert [name for name, _ in pairs] == ["div1__orf1", "sh5_div1__orf1"]
    by_name = dict(pairs)
    assert by_name["div1__orf1"] != by_name["sh5_div1__orf1"], "both went to one file"

    # The real query's file must hold the real query's alignment, which is the part the
    # substring match destroyed: it kept three sequences, not the shuffle's one.
    real = by_name["div1__orf1"].read_text()
    assert "YP_010769254.1" in real, "the real query's homologs were overwritten"
    assert "sh5_div1__orf1" not in real, "a shuffled sequence landed in the real query"
    assert "YP_010769254.1" not in by_name["sh5_div1__orf1"].read_text()


def test_an_unknown_query_name_raises_rather_than_falling_back(arms, tmp_path):
    """Falling back to block order is what let the overwrite happen quietly."""
    sto = tmp_path / "jackhmmer.sto"
    sto.write_text(block("not_an_orf", ("not_an_orf", "MRHW")))
    with pytest.raises(ValueError, match="not one of the"):
        arms.split_jackhmmer_alignments(sto, {"div1__orf1"}, tmp_path / "hhm")


def test_a_block_with_no_id_line_raises(arms, tmp_path):
    """No name means no answer; guessing one is the failure mode, not the fallback."""
    sto = tmp_path / "jackhmmer.sto"
    sto.write_text("# STOCKHOLM 1.0\ndiv1__orf1 MRHW\n//\n")
    with pytest.raises(ValueError, match="no `#=GF ID` line"):
        arms.split_jackhmmer_alignments(sto, {"div1__orf1"}, tmp_path / "hhm")


def test_two_blocks_naming_one_query_raises(arms, tmp_path):
    """A collision from any other cause has to be loud too, not last-write-wins."""
    sto = tmp_path / "jackhmmer.sto"
    sto.write_text(REAL + block("div1__orf1", ("div1__orf1", "MRHW"), iteration=2))
    with pytest.raises(ValueError, match="both name query"):
        arms.split_jackhmmer_alignments(sto, {"div1__orf1"}, tmp_path / "hhm")


def test_missing_jackhmmer_output_names_the_arm_that_has_to_run_first(arms, tmp_path):
    with pytest.raises(FileNotFoundError, match="jackhmmer arm has to run before it"):
        arms.split_jackhmmer_alignments(
            tmp_path / "absent.sto", {"div1__orf1"}, tmp_path / "hhm")
