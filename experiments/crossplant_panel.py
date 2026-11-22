"""Cross-plant subject panel: other wastewater plants' RNA libraries, plus a Berlin comparator.

Same columns as `results/recurrence/panel.tsv` (with project and plant appended), so
`recurrence_pass.py` can consume it. Built from NCBI runinfo and BioSample records; ENA has
no sample records for PRJNA1543313. Every filter is a fixed rule below, and a count that
does not match the expected one is a refusal.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import random
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from subtract_known import provenance  # noqa: E402

RUNINFO = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc={}"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
COLUMNS = ["run", "date", "window", "water", "prep", "capture_enriched", "paired",
           "read_count", "sra_key", "project", "plant"]
EXPECTED = {"PRJNA729801": 277, "PRJNA1543313": 26, "PRJNA924011": 55,
            "PRJNA842541": 140, "PRJNA941081": None}
QUERY_SAMPLE = "SRR23973193"
N_COMPARATOR = 6


def fetch(url: str, data: bytes | None = None, retries: int = 5) -> str:
    for i in range(retries):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "hapax-panel"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode()
        except Exception:
            if i == retries - 1:
                raise
            time.sleep(2 ** i)


def runinfo(project: str) -> list[dict]:
    rows = csv.DictReader(io.StringIO(fetch(RUNINFO.format(project))))
    return [r for r in rows if r.get("Run", "").startswith(("SRR", "ERR", "DRR"))]


def biosamples(accs: list[str]) -> dict[str, dict]:
    """BioSample accession -> {'title', 'name', attributes...}."""
    out = {}
    for i in range(0, len(accs), 200):
        body = f"db=biosample&id={','.join(accs[i:i + 200])}".encode()
        root = ET.fromstring(fetch(EFETCH, data=body))
        for s in root.iter("BioSample"):
            rec = {a.get("attribute_name"): (a.text or "") for a in s.iter("Attribute")}
            rec["_title"] = s.findtext("Description/Title") or ""
            out[s.get("accession")] = rec
        time.sleep(0.4)
    missing = set(accs) - set(out)
    if missing:
        sys.exit(f"refusing: {len(missing)} BioSamples not returned, e.g. {sorted(missing)[:3]}")
    return out


def collection_date(bs: dict) -> str:
    d = bs.get("collection_date", "").split("/")[0].strip()
    if len(d) != 10 or d[4] != "-":
        raise ValueError(f"unparsable collection_date {bs.get('collection_date')!r}")
    return d


def socal(r, bs):
    if r["LibraryStrategy"] != "RNA-Seq":
        return None
    plant = r["LibraryName"].split("_")[0]
    return plant, collection_date(bs), ""


def sixplants(r, bs):
    if "_RAW_" not in r["SampleName"]:
        return None
    f = r["SampleName"].split("_")
    return f[1], collection_date(bs), ""


def ohio(r, bs):
    if r["LibraryLayout"] != "PAIRED" or not r["Model"].startswith("NextSeq"):
        return None
    return r["SampleName"].split("_")[0], collection_date(bs), ""


def india(r, bs):
    return "", collection_date(bs), ""


def korea(r, bs):
    if r["LibraryStrategy"] != "RNA-Seq" or "influent" not in bs["_title"].lower():
        return None
    return bs["_title"].split()[0], collection_date(bs), ""


RULES = {"PRJNA729801": socal, "PRJNA1543313": sixplants, "PRJNA924011": ohio,
         "PRJNA842541": india, "PRJNA941081": korea}


def project_rows(project: str) -> list[dict]:
    runs = runinfo(project)
    bs = biosamples(sorted({r["BioSample"] for r in runs}))
    out = []
    for r in runs:
        got = RULES[project](r, bs[r["BioSample"]])
        if got is None:
            continue
        plant, date, window = got
        if project in ("PRJNA842541", "PRJNA941081"):  # several samples per date
            water = "_".join(x for x in (date, plant, r["SampleName"]) if x)
        else:
            water = "_".join(x for x in (date, plant, window) if x)
        out.append({"run": r["Run"], "date": date, "window": window, "water": water,
                    "prep": "", "capture_enriched": 0,
                    "paired": int(int(r["spots_with_mates"] or 0) > 0),
                    "read_count": int(r["spots"]), "project": project, "plant": plant,
                    "sra_key": f"s3://sra-pub-run-odp/sra/{r['Run']}/{r['Run']}"})
    want = EXPECTED[project]
    if want is not None and len(out) != want:
        sys.exit(f"refusing: {project} gave {len(out)} runs, expected {want}")
    if not out:
        sys.exit(f"refusing: {project} gave no runs")
    return out


def comparator(path: Path, seed: int) -> list[dict]:
    rows = list(csv.DictReader(path.open(), delimiter="\t"))
    own = next(r["water"] for r in rows if r["run"] == QUERY_SAMPLE)
    ok = [r for r in rows if r["prep"] == "" and r["capture_enriched"] == "0"
          and r["water"] != own]
    dates = sorted({r["date"] for r in ok})
    picks = [dates[round(i * (len(dates) - 1) / (N_COMPARATOR - 1))]
             for i in range(N_COMPARATOR)]
    assert len(set(picks)) == N_COMPARATOR
    rng = random.Random(seed)
    out = []
    for d in picks:
        r = rng.choice(sorted((x for x in ok if x["date"] == d), key=lambda x: x["run"]))
        out.append({**{c: r[c] for c in COLUMNS[:9]}, "project": "PRJNA948850",
                    "plant": "Berlin"})
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--berlin-panel", type=Path, default=Path("results/recurrence/panel.tsv"))
    ap.add_argument("--out", type=Path, default=Path("results/crossplant"))
    args = ap.parse_args()
    t0 = time.time()

    rows = []
    for p in RULES:
        rows += project_rows(p)
        time.sleep(0.5)
    rows += comparator(args.berlin_panel, args.seed)
    if len({r["run"] for r in rows}) != len(rows):
        sys.exit("refusing: duplicate run accessions")

    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out / "panel.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    per = {}
    for r in rows:
        p = per.setdefault(r["project"], {"runs": 0, "paired": 0, "plants": set(),
                                          "waters": set(), "dates": set(), "reads": 0})
        p["runs"] += 1
        p["paired"] += r["paired"] if isinstance(r["paired"], int) else int(r["paired"])
        p["plants"].add(r["plant"]); p["waters"].add(r["water"]); p["dates"].add(r["date"])
        p["reads"] += int(r["read_count"])
    summary = {k: {"runs": v["runs"], "paired": v["paired"], "plants": len(v["plants"]),
                   "waters": len(v["waters"]), "dates": len(v["dates"]),
                   "date_range": [min(v["dates"]), max(v["dates"])],
                   "reads_millions": round(v["reads"] / 1e6, 1)} for k, v in per.items()}
    manifest = {"script": "experiments/crossplant_panel.py", **provenance(__file__),
                "parameters": {"seed": args.seed, "berlin_panel": str(args.berlin_panel),
                               "n_comparator": N_COMPARATOR, "query_sample": QUERY_SAMPLE,
                               "dropped_projects": ["PRJNA902729", "PRJNA1445946"]},
                "source": "NCBI SRA runinfo + BioSample efetch",
                "fetched": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "seconds": round(time.time() - t0, 1),
                "rows": len(rows), "summary": summary}
    (args.out / "panel_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
