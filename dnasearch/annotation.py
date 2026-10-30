"""Pull the handful of annotation facts triage needs out of a Variant.

Field names differ between pipelines (nf-core/raredisease vcfanno+VEP vs. the
lightweight validation annotation), so each fact is looked up from a list of
"INFO:<key>" / "CSQ:<field>" sources configured in the case file.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .vcf import Variant

# Sequence Ontology consequence classes, most to least severe.
LOF = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "start_lost", "stop_lost",
    "transcript_amplification",
}
MODERATE = {
    "missense_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant", "splice_region_variant",
    "splice_donor_5th_base_variant", "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
}
LOW = {
    "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
    "intron_variant", "upstream_gene_variant", "downstream_gene_variant",
    "non_coding_transcript_exon_variant", "start_retained_variant",
    "stop_retained_variant", "incomplete_terminal_codon_variant",
    "coding_sequence_variant", "mature_miRNA_variant", "NMD_transcript_variant",
    "non_coding_transcript_variant", "intergenic_variant",
    "regulatory_region_variant", "TF_binding_site_variant",
}

CLINVAR_PATHOGENIC = ("pathogenic", "likely_pathogenic")

DEFAULT_SOURCES = {
    "af": ["INFO:gnomAD_AF_grpmax", "INFO:gnomAD_AF", "INFO:GNOMADAF_popmax",
           "INFO:GNOMADAF", "CSQ:gnomADg_AF", "CSQ:gnomAD_AF", "CSQ:MAX_AF"],
    "clinvar": ["INFO:CLNSIG", "INFO:CLNSIG_MOD", "CSQ:ClinVar_CLNSIG", "CSQ:CLIN_SIG"],
    # ClinVar GENEINFO ("SYMBOL:id|SYMBOL:id"): which gene(s) a ClinVar assertion is about
    "clinvar_gene": ["INFO:GENEINFO", "CSQ:ClinVar_GENEINFO"],
    "spliceai": ["CSQ:SpliceAI_pred_DS_max", "INFO:SpliceAI_DS_max"],
    "cadd": ["INFO:CADD", "CSQ:CADD_PHRED"],
}


@dataclass
class GeneEffect:
    gene: str
    consequence: str  # most severe term for this gene
    effect_class: str  # "lof" | "moderate" | "low"
    hgvsc: str = ""
    hgvsp: str = ""
    transcript: str = ""


@dataclass
class Annotated:
    variant: Variant
    effects: list[GeneEffect]
    af: float | None
    clinvar: str
    spliceai: float | None
    cadd: float | None
    clinvar_genes: set[str] = field(default_factory=set)

    def clinvar_applies_to(self, gene: str) -> bool:
        """ClinVar evidence counts for a gene only if ClinVar names it (when it names any)."""
        return not self.clinvar_genes or gene in self.clinvar_genes

    @property
    def clinvar_pathogenic(self) -> bool:
        c = self.clinvar.lower()
        # "Conflicting_classifications..." intentionally does not count
        return any(c.startswith(p) or f"/{p}" in c or f",{p}" in c
                   for p in CLINVAR_PATHOGENIC) and "conflicting" not in c

    @property
    def clinvar_benign(self) -> bool:
        c = self.clinvar.lower()
        return ("benign" in c and "pathogenic" not in c and "conflicting" not in c)


def _lookup(v: Variant, sources: list[str]) -> list[str]:
    vals = []
    for src in sources:
        kind, _, key = src.partition(":")
        if kind == "INFO" and key in v.info:
            vals.append(v.info[key])
        elif kind == "CSQ":
            vals.extend(c[key] for c in v.csq if c.get(key))
    return vals


def _max_float(vals: list[str]) -> float | None:
    best = None
    for raw in vals:
        for piece in raw.replace("&", ",").split(","):
            try:
                f = float(piece)
            except ValueError:
                continue
            best = f if best is None else max(best, f)
    return best


def effect_class(term: str) -> str:
    if term in LOF:
        return "lof"
    if term in MODERATE:
        return "moderate"
    return "low"


_RANK = {"lof": 0, "moderate": 1, "low": 2}


def annotate(v: Variant, sources: dict[str, list[str]] | None = None) -> Annotated:
    src = {**DEFAULT_SOURCES, **(sources or {})}
    per_gene: dict[str, GeneEffect] = {}
    for c in v.csq:
        gene = c.get("SYMBOL") or c.get("Gene") or ""
        if not gene:
            continue
        terms = (c.get("Consequence") or "").split("&")
        if not terms or terms == [""]:
            continue
        # most severe term in this transcript
        term = min(terms, key=lambda t: _RANK[effect_class(t)])
        eff = GeneEffect(gene, term, effect_class(term), c.get("HGVSc", ""),
                         c.get("HGVSp", ""), c.get("Feature", ""))
        prev = per_gene.get(gene)
        canonical = c.get("CANONICAL") == "YES" or bool(c.get("MANE_SELECT"))
        if (prev is None or _RANK[eff.effect_class] < _RANK[prev.effect_class]
                or (_RANK[eff.effect_class] == _RANK[prev.effect_class] and canonical)):
            per_gene[gene] = eff
    clin = _lookup(v, src["clinvar"])
    cv_genes = {g.split(":")[0] for raw in _lookup(v, src["clinvar_gene"])
                for g in raw.replace("&", "|").split("|") if g}
    return Annotated(
        variant=v,
        effects=list(per_gene.values()),
        af=_max_float(_lookup(v, src["af"])),
        clinvar=clin[0] if clin else "",
        spliceai=_max_float(_lookup(v, src["spliceai"])),
        cadd=_max_float(_lookup(v, src["cadd"])),
        clinvar_genes=cv_genes,
    )
