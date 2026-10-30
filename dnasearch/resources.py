"""Download / refresh reference knowledge (HPO, PanelApp panels)."""

from __future__ import annotations

import json
import urllib.request
from pathlib import Path

HPO_BASE = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/"
HPO_FILES = ["hp.obo", "genes_to_phenotype.txt", "genes_to_disease.txt", "phenotype.hpoa"]
PANELAPP = "https://panelapp.genomicsengland.co.uk/api/v1/panels/{id}/"


def _get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "dnasearch/0.1"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read()


def fetch_hpo(dest: str) -> None:
    d = Path(dest)
    d.mkdir(parents=True, exist_ok=True)
    for f in HPO_FILES:
        print(f"downloading {f}")
        (d / f).write_bytes(_get(HPO_BASE + f))


def fetch_panelapp(panel_id: int, out: str) -> None:
    data = json.loads(_get(PANELAPP.format(id=panel_id)))
    rows = []
    for g in data["genes"]:
        gd = g["gene_data"]
        loc = ""
        for build, recs in (gd.get("ensembl_genes") or {}).items():
            if build.lower() == "grch38":
                loc = next(iter(recs.values()))["location"]
        rows.append((gd["hgnc_symbol"], g["confidence_level"], g["mode_of_inheritance"],
                     gd.get("hgnc_id", ""), loc))
    rows.sort()
    with open(out, "w") as fh:
        fh.write(f"## source: Genomics England PanelApp panel {panel_id} '{data['name']}' "
                 f"version {data['version']} ({data.get('version_created', '')[:10]})\n")
        fh.write("## confidence: 3=green (diagnostic grade), 2=amber, 1=red\n")
        fh.write("symbol\tconfidence\tmode_of_inheritance\thgnc_id\tgrch38_location\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {len(rows)} genes to {out}")


def panel_bed(panel_tsv: str, out: str, padding: int = 5000, min_confidence: int = 3) -> None:
    """BED (chr-prefixed, GRCh38) of panel gene bodies plus padding."""
    from .genes import read_tsv
    recs = []
    for r in read_tsv(panel_tsv):
        if int(r["confidence"]) < min_confidence or not r.get("grch38_location"):
            continue
        chrom, _, span = r["grch38_location"].partition(":")
        s, e = (int(x) for x in span.split("-"))
        recs.append((f"chr{chrom}", max(0, s - 1 - padding), e + padding, r["symbol"]))
    order = {f"chr{c}": i for i, c in enumerate([*map(str, range(1, 23)), "X", "Y", "MT"])}
    recs.sort(key=lambda x: (order.get(x[0], 99), x[1]))
    with open(out, "w") as fh:
        for r in recs:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {len(recs)} regions to {out}")
