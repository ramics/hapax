"""Fetch the real within-host haplotype panels and background genomes for the
multi-swarm metagenome (`experiments/swarm_metagenome.py`).

Why real panels: a swarm simulated by our own coalescent (python/swarm.py)
encodes our own assumptions about where variation falls. Single-genome
amplification (SGA/SGS) and clonal sets from one host are the closest thing to
ground-truth haplotypes that exist: each sequence is one template from one
host at one time. Their variation (V-loop indels in HIV env, superinfection in
HCV, overlapping ORFs across all of HBV) is what the assembler has to face.

Panels (one swarm = one host, deduplicated, since identical sequences cannot be
told apart by any assembler and would make recovery impossible by construction):

| swarm | virus / region | source | why |
|---|---|---|---|
| hiv_env_ZM247F | HIV-1 C env, ~2.6-3.2 kb | EU166xxx, Zambia SGA | high diversity (~5%), V-loop length variation |
| hiv_env_ZM246F | HIV-1 C env | same study | second high-diversity env swarm, same subtype: sister-population chimaera risk |
| hiv_gagpol_2139 | HIV-1 5' half (gag-pol-vif-vpr), ~5.35 kb | KC156179-KC156209 | gag/pol and pol/vif dual-coding overlaps |
| hcv_Pt3P, hcv_Pt3S | HCV 5' half (core..NS3 incl. ARF/F), ~4.9 kb | MT416885-MT416975 | one host, two strains (primary and superinfecting), kept as two swarms: a join between them is a chimaera |
| hcv_Pt7 | HCV 5' half | MT417081-MT417136 | the most diverse single-strain HCV host in the study |
| hbv_p03 | HBV complete genome, 3215 nt | KY881820-KY881843 | every base in two ORFs or more; tail to 13% |
| hbv_p10 | HBV complete genome | KY881980-KY882003 | second HBV host: sister-population chimaera risk |

Background: E. coli K-12, Bacteroides fragilis, Staphylococcus epidermidis,
phages lambda, T4, crAssphage, phiX174, plus HIV-1 HXB2 as a single known
subtype-B genome sitting beside the subtype-C swarms.

Writes data/swarm_metagenome/{panels,background}/*.fa and appends to
data/MANIFEST.md.
"""

from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

PANELS = {
    "hiv_env_ZM247F": ('Human immunodeficiency virus 1[orgn] AND ZM247F AND env AND 2400:3300[slen]', r"isolate ZM247F"),
    "hiv_env_ZM246F": ('Human immunodeficiency virus 1[orgn] AND ZM246F AND env AND 2400:3300[slen]', r"isolate ZM246F"),
    "hiv_gagpol_2139": ("KC156179:KC156209[accn]", r"isolate 2139"),
    "hcv_Pt3P": ("MT416885:MT416975[accn]", r"isolate Pt3\.Primary"),
    "hcv_Pt3S": ("MT416885:MT416975[accn]", r"isolate Pt3\.Super"),
    "hcv_Pt7": ("MT417081:MT417136[accn]", r"isolate Pt7\."),
    "hbv_p03": ("KY881820:KY881843[accn]", r"isolate patient-03 "),
    "hbv_p10": ("KY881980:KY882003[accn]", r"isolate patient-10 "),
}
BACKGROUND = ["NC_000913.3", "NC_003228.3", "NC_004461.1", "NC_001416.1", "NC_000866.4", "NC_024711.1", "NC_001422.1", "K03455.1"]


def get(url: str) -> str:
    for i in range(5):
        try:
            return urllib.request.urlopen(url, timeout=120).read().decode()
        except Exception:
            time.sleep(2 + 3 * i)
    raise RuntimeError(url)


def search(term: str) -> list[str]:
    r = json.loads(get(f"{E}/esearch.fcgi?db=nuccore&retmode=json&retmax=1000&term=" + urllib.parse.quote(term)))
    return r["esearchresult"]["idlist"]


def fetch_fasta(ids: list[str]) -> list[tuple[str, str]]:
    out = []
    for i in range(0, len(ids), 200):
        txt = get(f"{E}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id=" + ",".join(ids[i:i + 200]))
        name, seq = None, []
        for line in txt.splitlines():
            if line.startswith(">"):
                if name:
                    out.append((name, "".join(seq).upper()))
                name, seq = line[1:], []
            elif line.strip():
                seq.append(line.strip())
        if name:
            out.append((name, "".join(seq).upper()))
        time.sleep(0.4)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "data" / "swarm_metagenome"))
    a = ap.parse_args()
    out = Path(a.out)
    (out / "panels").mkdir(parents=True, exist_ok=True)
    (out / "background").mkdir(parents=True, exist_ok=True)
    rows = []
    for swarm, (term, pat) in PANELS.items():
        recs = [(n, s) for n, s in fetch_fasta(search(term)) if re.search(pat, n)]
        uniq: dict[str, str] = {}
        for n, s in recs:
            if set(s) <= set("ACGT"):
                uniq.setdefault(s, n.split()[0])
        with open(out / "panels" / f"{swarm}.fa", "w") as f:
            for s, n in sorted(uniq.items(), key=lambda x: x[1]):
                f.write(f">{swarm}|{n}\n{s}\n")
        accs = sorted(n.split()[0] for n, _ in recs)
        rows.append(f"| swarm panel {swarm} | `{term}` ({accs[0]}..{accs[-1]}) | {time.strftime('%Y-%m-%d')} | {len(recs)} records, {len(uniq)} unique ACGT-only haplotypes |")
        print(swarm, len(recs), len(uniq))
    for acc in BACKGROUND:
        (n, s), = fetch_fasta([acc])
        (out / "background" / f"{acc}.fa").write_text(f">{acc}\n{s}\n")
        rows.append(f"| background {acc} | {n[len(acc) + 1:][:60]} | {time.strftime('%Y-%m-%d')} | {len(s)} nt |")
        print(acc, len(s))
    man = ROOT / "data" / "MANIFEST.md"
    txt = man.read_text()
    marker = "\n## Multi-swarm metagenome (experiments/fetch_swarm_panels.py)\n"
    if marker not in txt:
        txt += marker + "\n| Purpose | Accession / query | Downloaded | Notes |\n|---|---|---|---|\n" + "\n".join(rows) + "\n"
        man.write_text(txt)


if __name__ == "__main__":
    main()
