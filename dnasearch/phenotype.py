"""HPO phenotype-to-gene similarity (Resnik, best-match average).

For each patient term p, the best match against a gene is the information
content of the most informative ancestor of p that is also an ancestor of any
of the gene's annotated terms. The gene score is the mean over patient terms,
normalised by the patients' own term IC so it falls in 0..1.
"""

from __future__ import annotations

import math
from collections import defaultdict

from .genes import read_tsv

ROOT = "HP:0000001"
PHENOTYPIC_ABNORMALITY = "HP:0000118"


class Ontology:
    def __init__(self, obo_path: str):
        self.parents: dict[str, set[str]] = defaultdict(set)
        self.names: dict[str, str] = {}
        self.alt: dict[str, str] = {}
        self._anc: dict[str, frozenset[str]] = {}
        cur: dict | None = None
        with open(obo_path) as fh:
            for line in fh:
                line = line.strip()
                if line == "[Term]":
                    cur = {"is_a": []}
                elif line.startswith("[") and line.endswith("]"):
                    cur = None
                elif cur is not None and line == "" and "id" in cur:
                    self._add(cur)
                    cur = None
                elif cur is not None and ":" in line:
                    k, _, val = line.partition(": ")
                    if k == "is_a":
                        cur["is_a"].append(val.split(" ")[0])
                    elif k == "alt_id":
                        cur.setdefault("alt_id", []).append(val)
                    elif k in ("id", "name", "is_obsolete", "replaced_by"):
                        cur[k] = val
            if cur and "id" in cur:
                self._add(cur)

    def _add(self, t: dict) -> None:
        tid = t["id"]
        self.names[tid] = t.get("name", "")
        if t.get("is_obsolete") == "true":
            if t.get("replaced_by"):
                self.alt[tid] = t["replaced_by"]
            return
        self.parents[tid].update(t["is_a"])
        for a in t.get("alt_id", []):
            self.alt[a] = tid

    def canonical(self, tid: str) -> str:
        return self.alt.get(tid, tid)

    def ancestors(self, tid: str) -> frozenset[str]:
        tid = self.canonical(tid)
        if tid in self._anc:
            return self._anc[tid]
        out = {tid}
        stack = [tid]
        while stack:
            for p in self.parents.get(stack.pop(), ()):
                if p not in out:
                    out.add(p)
                    stack.append(p)
        self._anc[tid] = frozenset(out)
        return self._anc[tid]


class GeneMatcher:
    def __init__(self, ont: Ontology, genes_to_phenotype: str):
        self.ont = ont
        self.gene_closure: dict[str, set[str]] = defaultdict(set)
        for row in read_tsv(genes_to_phenotype):
            sym, hpo = row.get("gene_symbol"), row.get("hpo_id")
            if sym and hpo:
                self.gene_closure[sym] |= ont.ancestors(hpo)
        n = len(self.gene_closure)
        counts: dict[str, int] = defaultdict(int)
        for terms in self.gene_closure.values():
            for t in terms:
                counts[t] += 1
        self.ic = {t: -math.log(c / n) for t, c in counts.items()} if n else {}

    def known_genes(self) -> set[str]:
        return set(self.gene_closure)

    def score(self, patient_terms: list[str], gene: str) -> float:
        closure = self.gene_closure.get(gene)
        if not closure or not patient_terms:
            return 0.0
        num = den = 0.0
        for p in patient_terms:
            anc = self.ont.ancestors(p)
            den += max((self.ic.get(t, 0.0) for t in anc), default=0.0)
            num += max((self.ic.get(t, 0.0) for t in anc & closure), default=0.0)
        return num / den if den else 0.0

    def matched_terms(self, patient_terms: list[str], gene: str) -> list[str]:
        """Patient terms the gene's annotations cover exactly or more specifically."""
        closure = self.gene_closure.get(gene, set())
        return [self.ont.names.get(self.ont.canonical(p), p) for p in patient_terms
                if self.ont.canonical(p) in closure]
