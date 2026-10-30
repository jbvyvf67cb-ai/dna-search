"""Write triage results as Markdown (for people), TSV and JSON (for tools)."""

from __future__ import annotations

import csv
import json
from datetime import date
from pathlib import Path

from .annotation import Annotated
from .triage import Case, Finding, Result

TIER_TEXT = {
    "1": "Tier 1: panel gene, inheritance fits, strong variant evidence",
    "2": "Tier 2: weaker evidence in a panel gene, or strong evidence in another disease gene",
    "3": "Tier 3: research candidates (gene with no established disease, or weak evidence with a good phenotype match)",
}

UNCOVERED_CLASSES = [
    ("Structural variants / CNVs", "sv",
     "Run Manta/CNVnator/GATK-gCNV (nf-core/raredisease does this) and pass the SV VCF separately."),
    ("Repeat expansions", None, "ExpansionHunter/STRipy output must be reviewed separately."),
    ("Mosaic / somatic variants (e.g. UBA1/VEXAS, clonal haematopoiesis)", None,
     "Germline callers miss low-fraction variants; needs a somatic-mode caller (Mutect2) and ideally deeper sequencing."),
    ("Mitochondrial heteroplasmy", None, "Use the mitochondrial outputs of nf-core/raredisease."),
    ("Methylation / imprinting", None, "Needs long-read or methylation-array data."),
    ("Gene expression / splicing effects", None, "Needs RNA-seq (e.g. DROP: OUTRIDER/FRASER)."),
]


def _variant_cells(a: Annotated, gene: str, case: Case) -> dict:
    v = a.variant
    eff = next((e for e in a.effects if e.gene == gene), None)
    gts = {s: v.calls[s].raw_gt for s in case.ped.people}
    return {
        "variant": v.key,
        "consequence": eff.consequence if eff else "",
        "hgvsc": eff.hgvsc if eff else "",
        "hgvsp": eff.hgvsp if eff else "",
        "gnomad_af": "" if a.af is None else f"{a.af:.2e}",
        "clinvar": a.clinvar,
        "spliceai": "" if a.spliceai is None else f"{a.spliceai:.2f}",
        "genotypes": " ".join(f"{s}={g}" for s, g in gts.items()),
    }


def _md_table(rows: list[dict], cols: list[str]) -> str:
    if not rows:
        return "_None._\n"
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in rows:
        out.append("| " + " | ".join(str(r.get(c, "")).replace("|", "\\|") for c in cols) + " |")
    return "\n".join(out) + "\n"


def _finding_rows(fs: list[Finding], case: Case) -> list[dict]:
    rows = []
    for fd in fs:
        for i, a in enumerate(fd.variants):
            r = {"_tier": fd.tier, "tier": fd.tier, "gene": fd.gene, "model": fd.model, "evidence": fd.strength,
                 "gene_source": fd.gene_source, "phenotype_score": f"{fd.phenotype_score:.2f}",
                 "matched_hpo": "; ".join(fd.matched_terms), "note": fd.note,
                 **_variant_cells(a, fd.gene, case)}
            if i:  # second variant of a compound het: don't repeat gene-level cells
                for k in ("tier", "model", "evidence", "gene_source", "phenotype_score", "matched_hpo", "note"):
                    r[k] = "↳"
            rows.append(r)
    return rows


def write(case: Case, res: Result, outdir: str) -> dict[str, str]:
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    all_rows = _finding_rows(res.findings, case)
    cols = ["tier", "gene", "model", "evidence", "phenotype_score", "variant", "consequence",
            "hgvsp", "hgvsc", "gnomad_af", "clinvar", "spliceai", "genotypes", "gene_source",
            "matched_hpo", "note"]

    tsv = out / f"{case.id}.findings.tsv"
    with open(tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(r for r in all_rows)

    js = out / f"{case.id}.json"
    js.write_text(json.dumps({
        "case": case.id, "proband": case.ped.proband, "hpo": case.hpo,
        "funnel": res.funnel, "qc": res.qc, "warnings": case.warnings,
        "findings": all_rows,
        "secondary_findings": _finding_rows(res.secondary, case) if case.secondary_findings else None,
        "carriers": _finding_rows(res.carriers, case) if case.carrier_status else None,
    }, indent=2))

    md = [f"# Variant triage report: {case.id}", "",
          f"Generated {date.today().isoformat()} by dnasearch. **Research use only. Not a clinical "
          "result.** Any finding must be confirmed in a CLIA-certified lab and interpreted by a "
          "clinical geneticist or immunologist before it informs care.", ""]
    md += ["## Case", "",
           f"- Proband: `{case.ped.proband}` ({res.qc['declared_sex']}); "
           f"{'trio' if case.ped.is_trio else 'singleton/partial family'} analysis",
           f"- Panels: {', '.join(case.panel_names) or 'none'}",
           f"- Phenotype (HPO): {', '.join(case.hpo) or 'none given'}"
           + ("" if case.matcher else " (phenotype ranking unavailable)"), ""]
    if case.warnings:
        md += ["## Warnings", ""] + [f"- {w}" for w in case.warnings] + [""]

    md += ["## Filter funnel (proband's variants)", "", "| step | variants |", "|---|---|"]
    md += [f"| {k} | {n:,} |" for k, n in res.funnel] + [""]

    q = res.qc
    md += ["## QC", ""]
    if q["mendelian_error_rate"] is not None:
        md.append(f"- Mendelian errors: {q['mendelian_errors']:,} / {q['mendelian_sites_checked']:,} "
                  f"sites ({q['mendelian_error_rate']:.3%})")
    if q["proband_chrX_het_fraction"] is not None:
        md.append(f"- Proband chrX (non-PAR) heterozygous fraction: {q['proband_chrX_het_fraction']:.2f} "
                  f"over {q['proband_chrX_nonPAR_alt_sites']:,} sites"
                  + (f", inferred sex: {q['inferred_sex_from_chrX']}" if "inferred_sex_from_chrX" in q else ""))
    md.append("")

    short = ["gene", "model", "evidence", "phenotype_score", "variant", "consequence", "hgvsp",
             "gnomad_af", "clinvar", "genotypes", "note"]
    for t in ("1", "2", "3"):
        rows = [r for r in all_rows if r["_tier"] == t]
        md += [f"## {TIER_TEXT[t]}", "", _md_table(rows, short)]

    md += ["## Secondary findings (ACMG SF genes)", ""]
    if case.secondary_findings:
        md.append(_md_table(_finding_rows(res.secondary, case),
                            ["gene", "variant", "consequence", "hgvsp", "clinvar", "genotypes"]))
    else:
        md += ["_Not evaluated: disabled in the case file (`report.secondary_findings`)._", ""]

    md += ["## Carrier status (recessive genes, one copy)", ""]
    if case.carrier_status:
        md.append(_md_table(_finding_rows(res.carriers, case),
                            ["gene", "variant", "consequence", "hgvsp", "clinvar", "genotypes"]))
    else:
        md += ["_Not evaluated: disabled in the case file (`report.carrier_status`)._", ""]

    md += ["## What this analysis could NOT assess", ""]
    hard = sorted(g for g in case.difficult if g in case.panels)
    if hard:
        md += ["Panel genes that 30x short-read WGS resolves poorly (a negative result here means little):", "",
               _md_table([{"gene": g, **case.difficult[g]} for g in hard],
                         ["gene", "issue", "recommended_method"])]
    if res.low_coverage_genes:
        md += [f"Panel genes with mean coverage below {case.min_gene_coverage:g}x:", "",
               ", ".join(f"{g} ({c:.1f}x)" for g, c in res.low_coverage_genes), ""]
    md += ["Variant classes not covered by this input:", ""]
    for name, kind, how in UNCOVERED_CLASSES:
        if kind == "sv" and res.sv_records_seen:
            continue
        md.append(f"- **{name}**: {how}")
    md.append("")
    md += ["## How to read this", "",
           "- Tiers are a triage order for human review, not classifications. Most healthy genomes "
           "produce some Tier 2/3 entries.",
           "- `compound_het_unphased` and `dominant_inherited` from an unaffected parent are "
           "automatically moved down one tier.",
           "- Re-run periodically: about 250 new gene-disease links are published each year.", ""]

    mdp = out / f"{case.id}.report.md"
    mdp.write_text("\n".join(md))
    return {"markdown": str(mdp), "tsv": str(tsv), "json": str(js)}

