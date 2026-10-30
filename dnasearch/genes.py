"""Gene-level knowledge: diagnostic panels, inheritance modes, secondary-findings
genes, and genes short-read 30x WGS cannot reliably assess."""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

AD, AR, XL = "AD", "AR", "XL"
ALL_MODES = frozenset({AD, AR, XL})


@dataclass
class GeneInfo:
    symbol: str
    modes: frozenset[str] = ALL_MODES
    panels: set[str] = field(default_factory=set)
    note: str = ""


def modes_from_panelapp(moi: str) -> frozenset[str]:
    m = moi.upper()
    if m.startswith("X-LINKED"):
        return frozenset({XL})
    if m.startswith("BOTH"):
        return frozenset({AD, AR})
    if m.startswith("BIALLELIC"):
        return frozenset({AR})
    if m.startswith("MONOALLELIC"):
        return frozenset({AD})
    return ALL_MODES


def read_tsv(path: str | Path) -> list[dict[str, str]]:
    with open(path, newline="") as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("##")), delimiter="\t"))


def load_panel(path: str | Path, name: str, min_confidence: int = 3) -> dict[str, GeneInfo]:
    """Panel TSV columns: symbol, confidence, mode_of_inheritance[, ...]."""
    out = {}
    for row in read_tsv(path):
        if int(row.get("confidence") or 3) < min_confidence:
            continue
        out[row["symbol"]] = GeneInfo(row["symbol"],
                                      modes_from_panelapp(row.get("mode_of_inheritance", "")),
                                      {name})
    return out


def load_symbol_list(path: str | Path, column: str = "symbol") -> dict[str, dict[str, str]]:
    return {r[column]: r for r in read_tsv(path)}


# HPO inheritance terms -> mode
HPO_MODES = {
    "HP:0000006": AD, "HP:0000007": AR, "HP:0001417": XL,
    "HP:0001419": XL, "HP:0001423": XL,
}


def load_hpo_gene_modes(genes_to_disease: str | Path, phenotype_hpoa: str | Path) -> dict[str, frozenset[str]]:
    """Genome-wide gene -> inheritance modes, derived from HPO disease annotations.
    Any gene appearing here is treated as an established disease gene."""
    disease_modes: dict[str, set[str]] = {}
    with open(phenotype_hpoa) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("database_id"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 10 and f[10] == "I" and f[3] in HPO_MODES:
                disease_modes.setdefault(f[0], set()).add(HPO_MODES[f[3]])
    gene_modes: dict[str, set[str]] = {}
    for row in read_tsv(genes_to_disease):
        sym, dis = row.get("gene_symbol"), row.get("disease_id")
        if not sym or not dis or row.get("association_type", "MENDELIAN") != "MENDELIAN":
            continue
        gene_modes.setdefault(sym, set()).update(disease_modes.get(dis, set()))
    return {g: frozenset(m) or ALL_MODES for g, m in gene_modes.items()}
