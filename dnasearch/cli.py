"""dnasearch command line."""

from __future__ import annotations

import argparse
import sys

from . import report, resources, triage


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="dnasearch", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("report", help="triage an annotated VCF and write a report")
    r.add_argument("--case", required=True, help="case TOML (pedigree, HPO terms, panels, options)")
    r.add_argument("--vcf", required=True, help="joint-called, VEP-annotated VCF(.gz)")
    r.add_argument("--out", required=True, help="output directory")

    f = sub.add_parser("fetch-resources", help="download HPO files and refresh a PanelApp panel")
    f.add_argument("--hpo-dir", default="resources/hpo")
    f.add_argument("--panel-id", type=int, default=398,
                   help="PanelApp panel id (398 = primary immunodeficiency)")
    f.add_argument("--panel-out", default="resources/panels/primary_immunodeficiency.tsv")
    f.add_argument("--skip-hpo", action="store_true")

    b = sub.add_parser("panel-bed", help="BED of panel gene bodies (GRCh38)")
    b.add_argument("--panel", required=True)
    b.add_argument("--out", required=True)
    b.add_argument("--padding", type=int, default=5000)

    a = ap.parse_args(argv)
    if a.cmd == "report":
        case = triage.load_case(a.case)
        res = triage.run(case, a.vcf)
        paths = report.write(case, res, a.out)
        for k, n in res.funnel:
            print(f"{k:45s} {n:>10,}")
        tiers = {t: sum(1 for x in res.findings if x.tier == t) for t in "123"}
        print(f"findings: tier1={tiers['1']} tier2={tiers['2']} tier3={tiers['3']}")
        for k, p in paths.items():
            print(f"{k}: {p}")
    elif a.cmd == "fetch-resources":
        if not a.skip_hpo:
            resources.fetch_hpo(a.hpo_dir)
        resources.fetch_panelapp(a.panel_id, a.panel_out)
    elif a.cmd == "panel-bed":
        resources.panel_bed(a.panel, a.out, a.padding)
    return 0


if __name__ == "__main__":
    sys.exit(main())
