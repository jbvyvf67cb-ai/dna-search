"""Summarise hap.py results for the trio as a Markdown table."""

import csv
import sys
from pathlib import Path

bench = Path(sys.argv[1])
print("| sample | type | truth variants | recall | precision | F1 |")
print("|---|---|---|---|---|---|")
for s in ("HG002", "HG003", "HG004"):
    f = bench / f"{s}.summary.csv"
    if not f.exists():
        continue
    for r in csv.DictReader(open(f)):
        if r["Filter"] != "PASS":
            continue
        print(f"| {s} | {r['Type']} | {r['TRUTH.TOTAL']} | {float(r['METRIC.Recall']):.4f} | "
              f"{float(r['METRIC.Precision']):.4f} | {float(r['METRIC.F1_Score']):.4f} |")
