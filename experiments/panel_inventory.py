"""The cross-sample panel: which runs of PRJNA948850 can serve as recurrence subjects.

WHY THIS IS A SCRIPT AND NOT A PARAGRAPH. The recurrence measurement's denominator is the
panel, so "present in 12 samples" means nothing until the panel is pinned, and two mistakes
are available in it that would both read as biology.

THE FIRST MISTAKE: COUNTING THE CAPTURE RUNS AS SAMPLES. The study's 215 RNA runs are 116
paired `RNA_WW_*` libraries and 99 single-end `xh_*` runs, and the GEO metadata for an
`xh_` sample records `enrichment: xHYB adv. agents` - Qiagen's xHYB adventitious-agent
hybrid-capture panel applied to the same RNA. So an `xh_` run is the same water through a
different instrument, one that deliberately enriches sequence related to a list of known
agents. Recurrence measured across a mixture of the two would be partly a statement about
what that panel targets, enriched where a residue read happens to resemble something on it
and depleted elsewhere, and the direction would differ read by read. They are kept as a
separate instrument with their own question ("does capture pull our residue up?"), and the
panel is the 116.

THE SECOND MISTAKE: TREATING RUNS AS INDEPENDENT WATER. Of the 116, several share a
collection date and differ only in preparation - `cer_trizol`, `cer_BIOM`, `centri_BIOM`,
`PCRinh_rem`, and the five 220404 capture-workflow variants. A molecule present in one is
expected in its sibling for a reason that has nothing to do with recurring in the sewer, so
occupancy over runs is inflated relative to occupancy over dates by exactly the amount the
duplicated dates contribute. Both are reported; neither is corrected into the other.

WHAT THIS WRITES. One TSV row per RNA run with the fields a subject needs - accession, the
sample title it comes from, parsed date and collection window, preparation suffix, whether
it is paired, whether it is capture-enriched, read and base counts, the compressed bytes, and
the `sra-pub-run-odp` key - plus a manifest. The sizes matter because they are what the
in-region staging is priced from: the panel is 130 GB of fastq.gz, which is why this cannot
run in the container that writes this file.

THE DATE PARSE IS ASSERTED, NOT TRUSTED. Titles are `RNA_WW_<yymmdd>_<window>[_<prep>]`, and
a title that does not parse is a run whose date would silently become `unknown` and whose
occupancy-by-date column would then be wrong for a reason no reader could see. So an
unparsable title is a refusal.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from subtract_known import provenance

ROOT = Path(__file__).resolve().parents[1]

PROJECT = "PRJNA948850"
ENA = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession={acc}&result=read_run"
       "&fields=run_accession,sample_alias,sample_title,library_strategy,library_source,"
       "read_count,base_count,fastq_bytes,fastq_ftp,collection_date,first_public"
       "&format=tsv")

#: SRA's own normalized copy, us-east-1, one object per run. The recurrence pass stages from
#: here rather than from ENA because ENA serves this container at about 1 MB/s.
SRA_BUCKET = "sra-pub-run-odp"

TITLE = re.compile(r"^(?P<kind>RNA_WW|xh)_(?P<date>\d{6})_(?P<rest>.+)$")

#: The collection window, taken from the part of the title after the date. Everything left
#: over is the preparation suffix, which is what makes two runs of one date siblings.
WINDOWS = ("20-22hmix", "22-24hmix", "20-22h_mix", "22-24h_mix", "20-22h", "22-24h",
           "24hmix", "24h_mix", "12-14h", "16-18h", "22-24", "morning",
           "KW1", "KW2", "KW3", "KW6")

#: Preparation suffixes, and the list is asserted rather than inferred. Whatever is left of
#: a title after date and window is either a preparation variant - in which case two runs
#: sharing a date and window are the SAME WATER and must not count as two samples - or a
#: further piece of window information, like the `_1`/`_2` and `_22-24` that follow some
#: windows. Getting that split wrong changes the denominator of every occupancy figure, so an
#: unrecognised suffix is a refusal rather than a guess in either direction.
PREPS = ("cer_trizol", "cer_BIOM", "cer_BIOMcombi", "centri_BIOM", "centri_BIOMcombi",
         "PCRinh_rem", "cer_trizol_PCRinh_rem", "befo_captu", "captu_elu",
         "captu_elu_conc", "captu_ft",
         # the xh_ titles' own spellings of the same things
         "CeresTrizolElu", "CeresBiomRNA", "CeresBiomCombiRNA", "concBiomCombiRNA",
         "CeresTrizolElu_PCRinhibRemov", "PCRinhibRemov", "beforeCaptu", "CaptuElu",
         "CaptuEluConc", "CaptuFT", "mix_CeresTrizolElu")

#: Suffixes that are part of the collection window, not the preparation: the two numbered
#: 24-hour composites of 220623, and the morning/evening split of 210325's KW2.
WINDOW_TAILS = ("1", "2", "22-24", "morning")


def fetch(url: str, retries: int = 5) -> str:
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=120) as h:
                return h.read().decode()
        except Exception as exc:
            if attempt == retries - 1:
                raise SystemExit(f"ENA did not answer after {retries} tries: {exc}")
            time.sleep(5 * (attempt + 1))
    raise AssertionError("unreachable")


def parse_title(title: str) -> dict:
    """Split a sample title into kind, date, window and preparation suffix."""
    m = TITLE.match(title)
    if not m:
        raise SystemExit(
            f"cannot parse sample title {title!r}. The occupancy-by-date column is built "
            f"from this parse, so an unparsed title would quietly become a sample with no "
            f"date rather than an error. Add its shape to TITLE/WINDOWS.")
    rest = m.group("rest")
    window, tail = rest, ""
    for cand in sorted(WINDOWS, key=len, reverse=True):
        if rest == cand:
            window, tail = cand, ""
            break
        if rest.startswith(cand + "_"):
            window, tail = cand, rest[len(cand) + 1:]
            break
    else:
        raise SystemExit(f"no known collection window in {title!r}; add it to WINDOWS")
    prep = ""
    if tail:
        head = tail.split("_")[0]
        if tail in PREPS:
            prep = tail
        elif head in WINDOW_TAILS:
            window = f"{window}_{head}"
            prep = tail[len(head) + 1:]
            if prep and prep not in PREPS:
                raise SystemExit(f"unrecognised preparation {prep!r} in {title!r}; see PREPS")
        else:
            raise SystemExit(
                f"unrecognised suffix {tail!r} in {title!r}. It is either a preparation "
                f"variant (same water, so not a second sample) or part of the collection "
                f"window (different water), and the occupancy denominator differs by which. "
                f"Add it to PREPS or WINDOW_TAILS rather than letting it default.")
    yy, mm, dd = m.group("date")[:2], m.group("date")[2:4], m.group("date")[4:6]
    # The window as a comparable key: `20-22hmix` and `20-22h_mix` are the same window
    # spelled two ways, and the capture runs use the second spelling throughout.
    return {"kind": m.group("kind"), "date": f"20{yy}-{mm}-{dd}",
            "window": window, "window_key": window.replace("_", ""), "prep": prep}


def rows_from(tsv: str) -> list[dict]:
    out = []
    for r in csv.DictReader(tsv.splitlines(), delimiter="\t"):
        if r["library_strategy"] != "RNA-Seq":
            continue                                   # the 24 DNA runs are another thread
        meta = parse_title(r["sample_title"])
        paired = ";" in r["fastq_ftp"]
        out.append({
            "run": r["run_accession"],
            "gsm": r["sample_alias"],
            "title": r["sample_title"],
            "date": meta["date"],
            "window": meta["window"],
            "water": f"{meta['date']}_{meta['window_key']}",
            "prep": meta["prep"],
            "kind": meta["kind"],
            "capture_enriched": int(meta["kind"] == "xh"),
            "paired": int(paired),
            "read_count": int(r["read_count"]),
            "base_count": int(r["base_count"]),
            "mean_read_len": round(int(r["base_count"]) / int(r["read_count"]) /
                                   (2 if paired else 1)),
            "fastq_bytes": sum(int(b) for b in r["fastq_bytes"].split(";") if b),
            "sra_key": f"s3://{SRA_BUCKET}/sra/{r['run_accession']}/{r['run_accession']}",
        })
    return sorted(out, key=lambda x: (x["kind"], x["date"], x["window"], x["prep"]))


def summarise(rows: list[dict]) -> dict:
    panel = [r for r in rows if not r["capture_enriched"]]
    capture = [r for r in rows if r["capture_enriched"]]
    dates = Counter(r["date"] for r in panel)
    waters = Counter(r["water"] for r in panel)
    return {
        "runs_total": len(rows),
        "panel_runs": len(panel),
        "panel_paired": sum(r["paired"] for r in panel),
        # THREE DENOMINATORS, because they differ by 2.8x and the honest one depends on the
        # question. A run is one library. A `water` is one date and one two-hour composite
        # window, so runs sharing it are preparation variants of the same water and a
        # molecule in one is expected in its sibling for reasons that are not recurrence. A
        # date pools the evening's two windows, which are different composites of the same
        # day's sewage - independent samples of a correlated thing.
        "panel_waters": len(waters),
        "panel_waters_with_prep_variants": sum(1 for n in waters.values() if n > 1),
        "panel_runs_sharing_water": sum(n for n in waters.values() if n > 1),
        "panel_dates": len(dates),
        "panel_dates_with_siblings": sum(1 for d, n in dates.items() if n > 1),
        "panel_runs_on_shared_dates": sum(n for n in dates.values() if n > 1),
        "panel_pairs_millions": round(sum(r["read_count"] for r in panel) / 1e6, 1),
        "panel_fastq_gb": round(sum(r["fastq_bytes"] for r in panel) / 1e9, 1),
        "panel_mean_read_len": round(sum(r["mean_read_len"] for r in panel) / len(panel)),
        "capture_runs": len(capture),
        "capture_paired": sum(r["paired"] for r in capture),
        "capture_fastq_gb": round(sum(r["fastq_bytes"] for r in capture) / 1e9, 1),
        "prep_variants": dict(Counter(r["prep"] for r in panel if r["prep"])),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1,
                    help="kept for the repo's convention; this run is deterministic given "
                         "the archive's answer, and the answer itself is the seed")
    ap.add_argument("--project", default=PROJECT)
    ap.add_argument("--out", type=Path, default=ROOT / "results/recurrence")
    ap.add_argument("--query-sample", default="SRR23973193",
                    help="the sample the residue is defined on; it is excluded from its own "
                         "panel by the recurrence pass, and named here so the exclusion is "
                         "a recorded parameter rather than a line of code")
    args = ap.parse_args()

    t0 = time.perf_counter()
    rows = rows_from(fetch(ENA.format(acc=args.project)))
    if not any(r["run"] == args.query_sample for r in rows):
        raise SystemExit(f"{args.query_sample} is not in {args.project}'s RNA runs")
    args.out.mkdir(parents=True, exist_ok=True)

    fields = list(rows[0].keys())
    with open(args.out / "panel.tsv", "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    summary = summarise(rows)
    query = next(x for x in rows if x["run"] == args.query_sample)
    partner = [r["run"] for r in rows if r["capture_enriched"] and r["water"] == query["water"]]
    summary["query_water"] = query["water"]
    summary["query_prep_siblings"] = [r["run"] for r in rows if not r["capture_enriched"]
                                      and r["water"] == query["water"]
                                      and r["run"] != query["run"]]
    summary["query_sample"] = args.query_sample
    summary["query_capture_partner"] = partner
    manifest = {"script": "experiments/panel_inventory.py", **provenance(__file__),
                "parameters": {"project": args.project, "seed": args.seed,
                               "query_sample": args.query_sample},
                "source": "ENA filereport, read_run", "fetched": time.strftime(
                    "%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "seconds": round(time.perf_counter() - t0, 1), "summary": summary}
    (args.out / "panel_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
