"""Reading BLAST and DIAMOND tabular output, in one place, with its units stated.

WHY THIS FILE EXISTS. About fourteen scripts in `experiments/` each request their own
`--outfmt 6` column list, unpack it by position, and pick a best hit per query. They agree
on the core (strict `>` on bitscore, so the first-seen hit wins a tie) and drift on
everything around it:

- **Units.** For a translated search (blastx, DIAMOND blastx) `length` is in amino acids
  and `qlen`, `qstart` and `qend` are query nucleotides. `holdout_cascade.parse_best`
  divided one by the other, and every protein best-HSP coverage capped at exactly 0.3333.
  That looks like a plausible low coverage, not like a bug. `divergent_known` computed the
  same quantity as `3 * aln_aa / contig_nt`. Two copies differing by 3x is how a number
  becomes unreadable. Here the program is a required argument, and every length carries
  its unit in its name.
- **Absence.** A missing output file was read as "no hits" by more than one copy. That
  conflates "searched, zero hits" with "never searched", which this repo has paid for
  twice. Here a missing file raises and an empty one is an empty list.
- **Malformed rows.** Copies skipped a row with too few fields. A column list that
  disagrees with the command that wrote the file then reads as fewer hits. Here it raises.
- **The E-value.** Three copies dropped it. An identity without its E-value is not
  quotable under this repo's rule, so every row here carries it.

Nothing here filters by significance unless asked. A best hit is the best HSP by
bitscore, whatever its E-value; floors are the caller's operating point and belong in the
manifest beside the number.
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable, Iterable


class SearchOutputError(RuntimeError):
    """Raised instead of returning something that reads as a search result."""


# How each program's lengths are measured: (query unit, alignment-length unit, subject unit).
UNITS = {
    "blastn": ("nt", "nt", "nt"),
    "megablast": ("nt", "nt", "nt"),
    "blastp": ("aa", "aa", "aa"),
    "diamond_blastp": ("aa", "aa", "aa"),
    "blastx": ("nt", "aa", "aa"),
    "diamond_blastx": ("nt", "aa", "aa"),
    "tblastn": ("aa", "aa", "nt"),
}

_INT = {"length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "qlen", "slen",
        "nident", "positive", "gaps", "qframe", "sframe"}
_FLOAT = {"pident", "evalue", "bitscore", "score", "qcovhsp", "qcovs", "scovhsp", "ppos"}

# The column list most of this repo already requests; scripts that want more name them.
DEFAULT_COLUMNS = ("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                   "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen")


def outfmt_args(columns: Iterable[str] = DEFAULT_COLUMNS, *, diamond: bool) -> list[str]:
    """The `--outfmt`/`-outfmt` arguments that make a file `read_hits` can parse.

    The command and the parser take their column list from the same place, so the two
    cannot disagree about which position holds the E-value.
    """
    cols = list(columns)
    return ["--outfmt", "6", *cols] if diamond else ["-outfmt", " ".join(["6", *cols])]


def read_hits(path: Path, columns: Iterable[str] = DEFAULT_COLUMNS, *,
              program: str) -> list[dict]:
    """Every HSP in a tabular file, as a dict keyed by column name, with units attached.

    Raises if the file is missing (never searched is not zero hits), if any row has a
    field count different from `columns`, or if `program` is not one whose units are
    known. Every row also carries `program`, so a coverage computed later cannot forget
    which units it is in.
    """
    if program not in UNITS:
        raise SearchOutputError(f"unknown program {program!r}; units are known for "
                                f"{sorted(UNITS)}")
    cols = list(columns)
    for need in ("qseqid", "sseqid", "evalue", "bitscore"):
        if need not in cols:
            raise SearchOutputError(f"column list lacks {need!r}; a hit without it is "
                                    f"not interpretable")
    path = Path(path)
    if not path.exists():
        raise SearchOutputError(f"{path} does not exist: that is 'not searched', not "
                                f"'zero hits'")
    rows = []
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            if len(f) != len(cols):
                raise SearchOutputError(f"{path}:{ln} has {len(f)} fields, the column "
                                        f"list has {len(cols)}; the command and the "
                                        f"parser disagree")
            row = {"program": program}
            for c, v in zip(cols, f):
                row[c] = int(v) if c in _INT else float(v) if c in _FLOAT else v
            rows.append(row)
    return rows


def pair_key(qseqid: str) -> str:
    """`name/1` and `name/2` -> `name`: the read pair, this repo's counting unit."""
    return qseqid.rsplit("/", 1)[0]


def best_per_query(rows: Iterable[dict],
                   key: Callable[[str], str] | None = None) -> dict[str, dict]:
    """Best HSP per query by bitscore. On a tie the first-seen row wins.

    Best by bitscore and not by identity, because identity on a short alignment is
    quantised and rises as the alignment shortens. Ordering by it prefers the least
    informative hit, and this repo has the retraction to prove it. Strict `>` keeps the
    tie rule every existing copy has, so a script switching to this function picks the
    same rows. `key` maps a qseqid to the unit counted. Pass `pair_key` to count pairs.
    The row's own `qseqid` is kept, so which mate won is still recoverable.
    """
    best: dict[str, dict] = {}
    for r in rows:
        k = key(r["qseqid"]) if key else r["qseqid"]
        cur = best.get(k)
        if cur is None or r["bitscore"] > cur["bitscore"]:
            best[k] = r
    return best


def aligned_query_nt(row: dict) -> int:
    """Length of the HSP on the query, in query units, from the query coordinates.

    From `qstart`/`qend` rather than `length`, because `length` counts gaps and, under a
    translated search, counts residues. The coordinates are what the query actually
    spans. For blastx that is nucleotides, and so is `qlen`.
    """
    for c in ("qstart", "qend"):
        if c not in row:
            raise SearchOutputError(f"query coverage needs {c!r}; request it in the "
                                    f"column list rather than deriving it from `length`")
    return abs(row["qend"] - row["qstart"]) + 1


def query_coverage(row: dict) -> float:
    """Fraction of the query this one HSP spans, in the query's own units. Always <= 1."""
    if "qlen" not in row:
        raise SearchOutputError("query coverage needs 'qlen' from the search itself; a "
                                "length read from elsewhere is how an off-by-one in a "
                                "FASTA parser becomes a claim")
    if row["qlen"] <= 0:
        raise SearchOutputError(f"qlen {row['qlen']} for {row['qseqid']}")
    return aligned_query_nt(row) / row["qlen"]


def subject_coverage(row: dict) -> float:
    """Fraction of the subject this HSP spans, from subject coordinates when present.

    Falls back to `length / slen` only when the alignment and subject are in the same
    unit (true for every program in UNITS except tblastn) and `sstart`/`send` were not
    requested. That ratio counts gaps, so it can exceed 1 slightly; coordinates cannot.
    """
    if "slen" not in row:
        raise SearchOutputError("subject coverage needs 'slen'")
    if "sstart" in row and "send" in row:
        return (abs(row["send"] - row["sstart"]) + 1) / row["slen"]
    _, aln_u, sub_u = UNITS[row["program"]]
    if aln_u != sub_u:
        raise SearchOutputError(f"{row['program']}: alignment length is {aln_u} and the "
                                f"subject is {sub_u}; request sstart/send")
    return row["length"] / row["slen"]


def _merge(spans: list[tuple[int, int]]) -> list[list[int]]:
    merged: list[list[int]] = []
    for a, b in sorted(spans):
        if merged and a <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return merged


def union_query_coverage(rows: Iterable[dict], *, per_subject: bool,
                         evalue_below: float | None = None) -> dict[str, dict]:
    """Union of HSP intervals on each query, as a covered fraction of `qlen`.

    `per_subject=True` unions one subject at a time and keeps the subject that covers the
    most. That answers "does a single known protein explain this contig". `per_subject=
    False` unions over every subject. That answers "how much of it do known proteins
    explain between them". For a divergent multi-gene genome the two differ, because each
    gene drifts onto a different nearest relative. `evalue_below` is a floor applied
    before the union. Without one, a pile of weak chance HSPs inflates a union without
    limit, because each adds interval.
    """
    iv: dict[tuple, list] = {}
    qlen: dict[str, int] = {}
    subs: dict[str, set] = {}
    for r in rows:
        if evalue_below is not None and not r["evalue"] < evalue_below:
            continue
        q = r["qseqid"]
        aligned_query_nt(r)                                 # raises without coordinates
        a, b = sorted((r["qstart"], r["qend"]))
        iv.setdefault((q, r["sseqid"]) if per_subject else (q,), []).append((a, b))
        subs.setdefault(q, set()).add(r["sseqid"])
        if "qlen" not in r or r["qlen"] <= 0:
            raise SearchOutputError(f"union coverage needs a positive qlen for {q}")
        qlen[q] = r["qlen"]
    out: dict[str, dict] = {}
    for k, spans in iv.items():
        q = k[0]
        merged = _merge(spans)
        covered = sum(b - a + 1 for a, b in merged)
        row = {"covered_bases": covered, "covered_fraction": covered / qlen[q],
               "blocks": len(merged), "qlen": qlen[q]}
        if per_subject:
            row["subject"] = k[1]
            cur = out.get(q)
            if cur is None or row["covered_fraction"] > cur["covered_fraction"]:
                out[q] = row
        else:
            row["significant_subjects"] = len(subs[q])
            out[q] = row
    return out
