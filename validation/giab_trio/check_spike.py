"""Compare dnasearch output on the spiked trio against the expected outcomes."""

import csv
import json
import sys
from pathlib import Path


def main(truth_tsv: str, report_dir: str) -> int:
    js = json.loads(next(Path(report_dir).glob("*.json")).read_text())
    tiers: dict[str, set[tuple[str, str]]] = {}
    # continuation rows of a compound het carry "↳"; map them to the owning finding's model
    last_model = ""
    for r in js["findings"]:
        last_model = r["model"] if r["model"] != "↳" else last_model
        tiers.setdefault(r["variant"], set()).add((r["_tier"], last_model))
    sf = {r["variant"] for r in js.get("secondary_findings") or []}
    car = {r["variant"] for r in js.get("carriers") or []}

    ok_all = True
    print(f"{'scenario':20s} {'gene':6s} {'expected':10s} {'observed':40s} result")
    for t in csv.DictReader(open(truth_tsv), delimiter="\t"):
        v, exp = t["variant"], t["expected"]
        obs = sorted(tiers.get(v, set()))
        if exp == "tier1":
            ok = any(tr == "1" and m == t["model"] for tr, m in obs)
        elif exp == "not_tiered":
            ok = not any(tr in ("1", "2") for tr, _ in obs)
        elif exp == "carrier":
            ok = v in car and not any(tr == "1" for tr, _ in obs)
        elif exp == "secondary":
            ok = v in sf
        else:
            ok = False
        desc = ", ".join(f"T{tr}:{m}" for tr, m in obs) or "-"
        desc += " +SF" if v in sf else ""
        desc += " +carrier" if v in car else ""
        print(f"{t['scenario']:20s} {t['gene']:6s} {exp:10s} {desc:40s} {'PASS' if ok else 'FAIL'}")
        ok_all &= ok
    print("ALL SPIKE-IN CHECKS PASSED" if ok_all else "SPIKE-IN CHECKS FAILED")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
