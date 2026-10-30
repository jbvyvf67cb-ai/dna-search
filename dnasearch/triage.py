"""Filter funnel + tiering: genome-wide variants -> short ranked review list."""

from __future__ import annotations

import tomllib
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from . import inheritance as inh
from .annotation import Annotated, annotate
from .genes import (ALL_MODES, GeneInfo, load_hpo_gene_modes, load_panel,
                    load_symbol_list)
from .pedigree import Pedigree, from_config
from .phenotype import GeneMatcher, Ontology
from .vcf import HET, HOM_ALT, VcfReader

DEFAULT_FILTERS = {
    "min_gq": 20,
    "min_dp": 8,
    "max_af_dominant": 0.001,
    "max_af_recessive": 0.01,
    "max_af_clinvar_pathogenic": 0.1,
    "spliceai_min": 0.2,
    "phenotype_tier3_min": 0.5,
}

DOMINANT_MODELS = {inh.DE_NOVO, inh.AD_INHERITED, inh.AD_UNPHASED,
                   inh.XL_HEMI, inh.XL_FEMALE_HET}
DOWNGRADED_MODELS = {inh.COMPHET_UNPHASED}


@dataclass
class Case:
    id: str
    ped: Pedigree
    hpo: list[str]
    filters: dict
    annotation_sources: dict
    panels: dict[str, GeneInfo]
    panel_names: list[str]
    disease_gene_modes: dict[str, frozenset[str]]
    matcher: GeneMatcher | None
    acmg_sf: dict[str, dict[str, str]]
    difficult: dict[str, dict[str, str]]
    secondary_findings: bool
    carrier_status: bool
    coverage_bed: str | None = None
    min_gene_coverage: float = 15.0
    warnings: list[str] = field(default_factory=list)


def _resolve(p: str, base: Path) -> Path:
    cand = base / p
    return cand if cand.exists() else Path(p)


def load_case(path: str) -> Case:
    base = Path(path).resolve().parent
    cfg = tomllib.loads(Path(path).read_text())
    res = cfg.get("resources", {})
    rep = cfg.get("report", {})
    warnings: list[str] = []

    panels: dict[str, GeneInfo] = {}
    panel_names = []
    for p in res.get("panels", []):
        panel_names.append(p["name"])
        for sym, gi in load_panel(_resolve(p["path"], base), p["name"],
                                  p.get("min_confidence", 3)).items():
            if sym in panels:
                panels[sym].panels |= gi.panels
            else:
                panels[sym] = gi

    disease_modes: dict[str, frozenset[str]] = {}
    matcher = None
    hpo_dir = res.get("hpo_dir")
    if hpo_dir and _resolve(hpo_dir, base).exists():
        d = _resolve(hpo_dir, base)
        if (d / "genes_to_disease.txt").exists() and (d / "phenotype.hpoa").exists():
            disease_modes = load_hpo_gene_modes(d / "genes_to_disease.txt", d / "phenotype.hpoa")
        if (d / "hp.obo").exists() and (d / "genes_to_phenotype.txt").exists():
            matcher = GeneMatcher(Ontology(str(d / "hp.obo")), str(d / "genes_to_phenotype.txt"))
    else:
        warnings.append("HPO resources not found: no phenotype ranking and no genome-wide "
                        "disease-gene list (run `dnasearch fetch-resources`).")

    case = cfg.get("case", {})
    return Case(
        id=case.get("id", Path(path).stem),
        ped=from_config({"samples": cfg["samples"], "proband": case.get("proband")}),
        hpo=case.get("hpo", []),
        filters={**DEFAULT_FILTERS, **cfg.get("filters", {})},
        annotation_sources=cfg.get("annotation", {}),
        panels=panels,
        panel_names=panel_names,
        disease_gene_modes=disease_modes,
        matcher=matcher,
        acmg_sf=load_symbol_list(_resolve(res["acmg_sf"], base)) if res.get("acmg_sf") else {},
        difficult=load_symbol_list(_resolve(res["difficult_genes"], base)) if res.get("difficult_genes") else {},
        secondary_findings=bool(rep.get("secondary_findings", False)),
        carrier_status=bool(rep.get("carrier_status", False)),
        coverage_bed=str(_resolve(res["coverage_bed"], base)) if res.get("coverage_bed") else None,
        min_gene_coverage=float(rep.get("min_gene_coverage", 15)),
        warnings=warnings,
    )


@dataclass
class Finding:
    tier: str
    gene: str
    model: str
    variants: list[Annotated]
    strength: str
    phenotype_score: float
    matched_terms: list[str]
    gene_source: str
    note: str = ""


@dataclass
class Result:
    funnel: list[tuple[str, int]]
    qc: dict
    findings: list[Finding]
    secondary: list[Finding]
    carriers: list[Finding]
    sv_records_seen: int
    low_coverage_genes: list[tuple[str, float]]


def variant_strength(a: Annotated, gene: str, f: dict, model: str) -> str | None:
    """'strong' | 'moderate' | None (not interesting)."""
    if a.clinvar_benign:
        return None
    eff = next((e for e in a.effects if e.gene == gene), None)
    if eff is None or not inh.in_gene(eff):
        return None  # ClinVar/consequence evidence belongs to the gene the variant is in
    if a.clinvar_pathogenic and a.clinvar_applies_to(gene):
        return "strong"
    cls = eff.effect_class
    if cls == "lof":
        return "strong"
    if cls == "moderate" and model == inh.DE_NOVO:
        return "strong"
    if cls == "moderate" or (a.spliceai or 0) >= f["spliceai_min"]:
        return "moderate"
    return None


def _af_ok(a: Annotated, model: str, f: dict) -> bool:
    af = a.af or 0.0
    if a.clinvar_pathogenic:
        return af <= f["max_af_clinvar_pathogenic"]
    lim = f["max_af_dominant"] if model in DOMINANT_MODELS else f["max_af_recessive"]
    return af <= lim


def _tier_down(t: str) -> str:
    return {"1": "2", "2": "3", "3": "3"}[t]


def run(case: Case, vcf_path: str) -> Result:
    f = case.filters
    rdr = VcfReader(vcf_path)
    case.ped.validate(rdr.samples)
    ped = case.ped
    pb = ped.proband
    counts = Counter()
    me_checked = me_errors = 0
    x_alt = x_het = 0
    sv_seen = 0
    no_af = 0
    cands: list[Annotated] = []

    for v in rdr:
        if v.alt.startswith("<") or "SVTYPE" in v.info:
            sv_seen += 1
        call = v.calls[pb]
        # trio QC on autosomal high-quality sites
        if ped.is_trio and v.passed and inh.norm_chrom(v.chrom) not in ("X", "Y", "M", "MT"):
            trio = [call, v.calls[ped.father.id], v.calls[ped.mother.id]]
            if all((c.gq or 0) >= f["min_gq"] for c in trio) and any(c.has_alt for c in trio):
                err = inh.mendelian_error(*trio)
                if err is not None:
                    me_checked += 1
                    me_errors += err
        if not call.has_alt:
            continue
        counts["proband variants"] += 1
        if inh.in_nonpar_x(v) and v.passed:
            x_alt += 1
            x_het += call.gt == HET
        if not (v.passed and (call.gq or 0) >= f["min_gq"] and (call.dp or 0) >= f["min_dp"]):
            continue
        counts["pass quality"] += 1
        a = annotate(v, case.annotation_sources)
        if not (any(e.effect_class != "low" for e in a.effects) or a.clinvar_pathogenic
                or (a.spliceai or 0) >= f["spliceai_min"]):
            continue
        counts["protein/splice-affecting or ClinVar P/LP"] += 1
        no_af += a.af is None
        if a.clinvar_benign:
            continue
        loosest = f["max_af_clinvar_pathogenic"] if a.clinvar_pathogenic else f["max_af_recessive"]
        if (a.af or 0.0) > loosest:
            continue
        counts["rare (population frequency)"] += 1
        cands.append(a)
    rdr.close()

    hits = inh.find_hits(cands, ped)
    by_key = {a.variant.key: a for a in cands}
    findings: list[Finding] = []
    secondary: list[Finding] = []
    carriers: list[Finding] = []
    seen = set()
    inherit_ok = set()

    for key, hs in hits.items():
        a = by_key[key]
        for h in hs:
            vs = [a] + ([by_key[h.partner]] if h.partner else [])
            if not all(_af_ok(x, h.model, f) for x in vs):
                continue
            strengths = [variant_strength(x, h.gene, f, h.model) for x in vs]
            if None in strengths:
                continue
            strength = "strong" if all(s == "strong" for s in strengths) else "moderate"
            dedup = (h.gene, h.model, tuple(sorted(x.variant.key for x in vs)))
            if dedup in seen:
                continue
            seen.add(dedup)
            inherit_ok.add(key)

            if h.gene in case.panels:
                modes, source = case.panels[h.gene].modes, "panel"
            elif h.gene in case.disease_gene_modes:
                modes, source = case.disease_gene_modes[h.gene], "known disease gene"
            else:
                modes, source = ALL_MODES, "no known disease"
            mode_ok = any(h.model in inh.MODE_MODELS[m] for m in modes)
            pscore = case.matcher.score(case.hpo, h.gene) if case.matcher else 0.0
            matched = case.matcher.matched_terms(case.hpo, h.gene) if case.matcher else []

            tier = None
            if mode_ok and source == "panel":
                tier = "1" if strength == "strong" else "2"
            elif mode_ok and source == "known disease gene":
                if strength == "strong":
                    tier = "2"
                elif pscore >= f["phenotype_tier3_min"]:
                    tier = "3"
            elif source == "no known disease" and strength == "strong" and \
                    h.model in (inh.DE_NOVO, inh.AR_HOM, inh.COMPHET, inh.XL_HEMI):
                tier = "3"
            note = h.note
            if tier and (h.model in DOWNGRADED_MODELS or "reduced penetrance" in h.note):
                tier = _tier_down(tier)
            if tier:
                findings.append(Finding(tier, h.gene, h.model, vs, strength, pscore,
                                        matched, source, note))

    # secondary findings and carrier status are evaluated independently of tiering
    if case.secondary_findings:
        secondary = secondary_findings(case, cands)
    for a in cands:
        for e in a.effects:
            if case.carrier_status and inh.in_gene(e) and a.variant.calls[pb].gt == HET and \
                    ((a.clinvar_pathogenic and a.clinvar_applies_to(e.gene)) or e.effect_class == "lof"):
                modes = (case.panels[e.gene].modes if e.gene in case.panels
                         else case.disease_gene_modes.get(e.gene))
                if modes == frozenset({"AR"}) and not any(
                        x.gene == e.gene and x.model in (inh.COMPHET, inh.AR_HOM) for x in findings):
                    carriers.append(Finding("carrier", e.gene, "het carrier", [a], "strong",
                                            0.0, [], "AR gene"))

    counts["fits an inheritance model"] = len(inherit_ok)
    counts["reported (tiers 1-3)"] = len({x.variant.key for fd in findings for x in fd.variants})
    order = ["proband variants", "pass quality", "protein/splice-affecting or ClinVar P/LP",
             "rare (population frequency)", "fits an inheritance model", "reported (tiers 1-3)"]
    findings.sort(key=lambda x: (x.tier, -x.phenotype_score, x.strength != "strong", x.gene))

    qc = {
        "mendelian_sites_checked": me_checked,
        "mendelian_errors": me_errors,
        "mendelian_error_rate": (me_errors / me_checked) if me_checked else None,
        "proband_chrX_nonPAR_alt_sites": x_alt,
        "proband_chrX_het_fraction": (x_het / x_alt) if x_alt else None,
        "declared_sex": {1: "male", 2: "female"}.get(ped.proband_person.sex, "unknown"),
    }
    affecting = counts["protein/splice-affecting or ClinVar P/LP"]
    qc["affecting_without_population_af"] = no_af
    # most protein-affecting variants in any genome are common; if many lack a frequency the
    # annotation is incomplete and the rarity filter is silently letting everything through
    if affecting >= 50 and no_af / affecting > 0.3:
        case.warnings.append(
            f"{no_af} of {affecting} protein-affecting variants have no population allele frequency; "
            "check the gnomAD annotation (field names in [annotation].af) before trusting rarity filtering")
    if qc["proband_chrX_het_fraction"] is not None and x_alt >= 50:
        inferred = "male" if qc["proband_chrX_het_fraction"] < 0.2 else "female"
        qc["inferred_sex_from_chrX"] = inferred
        if qc["declared_sex"] != "unknown" and inferred != qc["declared_sex"]:
            case.warnings.append(f"declared sex {qc['declared_sex']} but chrX heterozygosity suggests {inferred}")

    return Result([(k, counts[k]) for k in order], qc, findings, secondary, carriers,
                  sv_seen, low_coverage(case))


def secondary_findings(case: Case, cands: list[Annotated]) -> list[Finding]:
    """ACMG SF genes: ClinVar P/LP only, honouring each gene's allele-count and special rules."""
    pb = case.ped.proband_person
    per_gene: dict[str, list[tuple[Annotated, int]]] = {}
    for a in cands:
        if not a.clinvar_pathogenic:
            continue
        gt = a.variant.calls[pb.id].gt
        n = 2 if gt == HOM_ALT or (pb.sex == 1 and inh.in_nonpar_x(a.variant)) else 1
        for e in a.effects:
            rule = case.acmg_sf.get(e.gene) if inh.in_gene(e) else None
            if not rule:
                continue
            special = rule.get("special", "")
            if not a.clinvar_applies_to(e.gene):
                continue
            if special == "lof_only" and e.effect_class != "lof":
                continue
            if special == "homozygous_only" and gt != HOM_ALT:
                continue
            per_gene.setdefault(e.gene, []).append((a, n))
    out = []
    for gene, items in per_gene.items():
        rule = case.acmg_sf[gene]
        if sum(n for _, n in items) < int(rule.get("min_alleles") or 1):
            continue
        note = rule.get("report_rule", "")
        if rule.get("special") == "needs_case_review":
            note += " (needs expert review before reporting)"
        out.append(Finding("SF", gene, rule.get("inheritance", ""), [a for a, _ in items],
                           "strong", 0.0, [], f"ACMG SF: {rule.get('category', '')}", note))
    return sorted(out, key=lambda f: f.gene)


def low_coverage(case: Case) -> list[tuple[str, float]]:
    """Panel genes below min mean coverage, from a mosdepth --by <bed with gene names> regions file."""
    if not case.coverage_bed:
        return []
    from .vcf import open_text
    sums: dict[str, list[float]] = {}
    with open_text(case.coverage_bed) as fh:
        for line in fh:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            gene, length, cov = c[3], int(c[2]) - int(c[1]), float(c[4])
            s = sums.setdefault(gene, [0.0, 0.0])
            s[0] += cov * length
            s[1] += length
    out = [(g, s[0] / s[1]) for g, s in sums.items()
           if s[1] and g in case.panels and s[0] / s[1] < case.min_gene_coverage]
    return sorted(out, key=lambda x: x[1])
