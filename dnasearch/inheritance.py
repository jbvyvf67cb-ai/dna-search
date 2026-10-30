"""Inheritance-model matching for singleton or trio cases."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from .annotation import Annotated
from .genes import AD, AR, XL
from .pedigree import Pedigree, Person
from .vcf import HET, HOM_ALT, HOM_REF, Call, Variant

# GRCh38 pseudoautosomal regions
PAR = {"X": [(10_001, 2_781_479), (155_701_383, 156_030_895)],
       "Y": [(10_001, 2_781_479), (56_887_903, 57_217_415)]}

DE_NOVO = "de_novo"
AD_INHERITED = "dominant_inherited"
AD_UNPHASED = "dominant"  # singleton / parents missing
AR_HOM = "recessive_homozygous"
COMPHET = "compound_het"
COMPHET_UNPHASED = "compound_het_unphased"
XL_HEMI = "x_linked_hemizygous"
XL_FEMALE_HET = "x_linked_female_het"
XL_FEMALE_HOM = "x_linked_female_hom"

# which models satisfy which gene inheritance mode
MODE_MODELS = {
    AD: {DE_NOVO, AD_INHERITED, AD_UNPHASED},
    AR: {AR_HOM, COMPHET, COMPHET_UNPHASED},
    XL: {XL_HEMI, XL_FEMALE_HET, XL_FEMALE_HOM, DE_NOVO},
}


def norm_chrom(c: str) -> str:
    return c[3:] if c.lower().startswith("chr") else c


def in_nonpar_x(v: Variant) -> bool:
    c = norm_chrom(v.chrom)
    return c == "X" and not any(a <= v.pos <= b for a, b in PAR["X"])


NOT_IN_GENE = {"upstream_gene_variant", "downstream_gene_variant", "intergenic_variant"}


def in_gene(eff) -> bool:
    return eff.consequence not in NOT_IN_GENE


@dataclass
class Hit:
    gene: str
    model: str
    partner: str = ""  # other variant key for compound hets
    note: str = ""


def _gt(v: Variant, p: Person | None) -> int | None:
    if p is None:
        return None
    c = v.calls.get(p.id)
    return c.gt if c else None


def origin(v: Variant, ped: Pedigree) -> str:
    """Where the proband's alt allele came from: paternal/maternal/de_novo/both/unknown."""
    fa, mo = _gt(v, ped.father), _gt(v, ped.mother)
    if fa is None or mo is None:
        return "unknown"
    f_alt, m_alt = fa in (HET, HOM_ALT), mo in (HET, HOM_ALT)
    if not f_alt and not m_alt:
        return "de_novo"
    if f_alt and m_alt:
        return "both"
    return "paternal" if f_alt else "maternal"


def variant_models(a: Annotated, ped: Pedigree) -> list[str]:
    """Single-variant models (compound hets are handled per gene)."""
    v = a.variant
    pb = ped.proband_person
    g = _gt(v, pb)
    if g not in (HET, HOM_ALT):
        return []
    ori = origin(v, ped)
    models: list[str] = []
    unaffected = lambda p: p is not None and p.affected is False  # noqa: E731

    if in_nonpar_x(v):
        if pb.sex == 1:
            # a carrier mother is expected; an unaffected father with the allele argues against it
            if not (unaffected(ped.father) and _gt(v, ped.father) in (HET, HOM_ALT)):
                models.append(XL_HEMI)
            if ped.is_trio and _gt(v, ped.mother) == HOM_REF:
                models.append(DE_NOVO)
            return models
        if pb.sex == 2:
            models.append(XL_FEMALE_HOM if g == HOM_ALT else XL_FEMALE_HET)
    if norm_chrom(v.chrom) in ("Y", "M", "MT"):
        return models

    if ori == "de_novo":
        models.append(DE_NOVO)
    # an unaffected parent homozygous for the allele rules out a dominant effect
    unaffected_hom = any(unaffected(p) and _gt(v, p) == HOM_ALT for p in (ped.father, ped.mother))
    if g == HET and not unaffected_hom:
        if ori in ("paternal", "maternal"):
            # allele from an unaffected parent needs reduced penetrance; flagged in find_hits
            models.append(AD_INHERITED)
        elif ori in ("unknown", "both"):
            models.append(AD_UNPHASED)
    if g == HOM_ALT:
        if not (unaffected(ped.father) and _gt(v, ped.father) == HOM_ALT) and \
           not (unaffected(ped.mother) and _gt(v, ped.mother) == HOM_ALT):
            models.append(AR_HOM)
    return models


def find_hits(cands: list[Annotated], ped: Pedigree) -> dict[str, list[Hit]]:
    """Return variant key -> inheritance hits (one per gene/model)."""
    hits: dict[str, list[Hit]] = defaultdict(list)
    by_gene: dict[str, list[Annotated]] = defaultdict(list)
    pb = ped.proband_person
    for a in cands:
        models = variant_models(a, ped)
        for eff in a.effects:
            if not in_gene(eff):
                continue  # up/downstream of a gene says nothing about that gene
            for m in models:
                note = ""
                if m == AD_INHERITED:
                    ori = origin(a.variant, ped)
                    parent = ped.father if ori == "paternal" else ped.mother
                    note = f"inherited from {ori} parent" + (
                        " (unaffected: needs reduced penetrance)" if parent and parent.affected is False else "")
                if m == AD_UNPHASED and origin(a.variant, ped) == "both" and any(
                        p.affected is False for p in (ped.father, ped.mother)):
                    note = "both parents carry it, at least one unaffected: needs reduced penetrance"
                if m == AR_HOM and ped.is_trio and HOM_REF in (_gt(a.variant, ped.father), _gt(a.variant, ped.mother)):
                    note = "a parent is hom-ref: possible deletion on other allele or uniparental disomy"
                if m == XL_HEMI and _gt(a.variant, pb) == HET:
                    note = "heterozygous call in a male on chrX: check for mosaicism or mapping artefact"
                hits[a.variant.key].append(Hit(eff.gene, m, note=note))
            # compound-het candidates: heterozygous, not male non-PAR X
            if _gt(a.variant, pb) == HET and not (pb.sex == 1 and in_nonpar_x(a.variant)) \
                    and norm_chrom(a.variant.chrom) not in ("Y", "M", "MT"):
                by_gene[eff.gene].append(a)

    for gene, vs in by_gene.items():
        if len(vs) < 2:
            continue
        for i, a in enumerate(vs):
            for b in vs[i + 1:]:
                oa, ob = origin(a.variant, ped), origin(b.variant, ped)
                if _parent_excludes_pair(a.variant, b.variant, ped):
                    continue
                if not ped.is_trio or "unknown" in (oa, ob):
                    model, note = COMPHET_UNPHASED, "phase unknown: could be on the same copy (cis)"
                elif oa == ob and oa in ("paternal", "maternal"):
                    continue  # both from one parent: cis, not a compound het
                elif "both" in (oa, ob) and oa == ob:
                    model, note = COMPHET_UNPHASED, "both parents carry both variants: phase unresolved"
                else:
                    model, note = COMPHET, f"{oa} + {ob}"
                hits[a.variant.key].append(Hit(gene, model, partner=b.variant.key, note=note))
                hits[b.variant.key].append(Hit(gene, model, partner=a.variant.key, note=note))
    return hits


def _parent_excludes_pair(a: Variant, b: Variant, ped: Pedigree) -> bool:
    """An unaffected parent who is homozygous for either variant, or carries both, rules out
    the pair: either the parent would be affected, or the two variants travel together (cis)."""
    for p in (ped.father, ped.mother):
        if p is None or p.affected is not False:
            continue
        ga, gb = _gt(a, p), _gt(b, p)
        if HOM_ALT in (ga, gb):
            return True
        if ga in (HET, HOM_ALT) and gb in (HET, HOM_ALT):
            return True
    return False


def mendelian_error(child: Call, father: Call, mother: Call, child_male_x: bool = False) -> bool | None:
    """True if the trio genotype is Mendelian-inconsistent (autosomal), None if uncallable."""
    if None in (child.gt, father.gt, mother.gt) or child_male_x:
        return None
    alleles = {HOM_REF: {0}, HET: {0, 1}, HOM_ALT: {1}}
    possible = {f + m for f in alleles[father.gt] for m in alleles[mother.gt]}
    return child.gt not in possible
