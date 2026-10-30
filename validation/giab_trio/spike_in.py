"""Insert known ClinVar pathogenic variants into the GIAB trio VCF with chosen
genotypes, so the triage logic can be checked against known answers.

Each scenario picks the first ClinVar variant in the gene that is a
multi-submitter/expert-panel 'Pathogenic' SNV with a missense or nonsense
consequence and is not already a site in the trio.
"""

import argparse
import gzip
import re

# (scenario, gene, genotypes HG002/HG003/HG004, expected outcome, expected model)
SCENARIOS = [
    ("x_linked", "BTK", [("1/1", "0/0", "0/1")], "tier1", "x_linked_hemizygous"),
    ("de_novo_dominant", "STAT3", [("0/1", "0/0", "0/0")], "tier1", "de_novo"),
    ("compound_het", "RAG1", [("0/1", "0/1", "0/0"), ("0/1", "0/0", "0/1")], "tier1", "compound_het"),
    # both variants from the father = same chromosome copy: must NOT be called compound het
    ("cis_negative", "RAG2", [("0/1", "0/1", "0/0"), ("0/1", "0/1", "0/0")], "not_tiered", ""),
    ("carrier_only", "DOCK8", [("0/1", "0/0", "0/1")], "carrier", ""),
    ("secondary_finding", "BRCA2", [("0/1", "0/1", "0/0")], "secondary", ""),
]
GOOD_REVIEW = ("criteria_provided,_multiple_submitters", "reviewed_by_expert_panel", "practice_guideline")


def clinvar_candidates(path):
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            ref, alt = f[3], f[4]
            if len(ref) != 1 or len(alt) != 1:
                continue
            info = f[7]
            gene = re.search(r"GENEINFO=([^:;]+)", info)
            sig = re.search(r"CLNSIG=([^;]+)", info)
            rev = re.search(r"CLNREVSTAT=([^;]+)", info)
            mc = re.search(r"MC=([^;]+)", info)
            if not (gene and sig and rev and mc) or sig.group(1) != "Pathogenic":
                continue
            if not any(r in rev.group(1) for r in GOOD_REVIEW):
                continue
            if "missense_variant" not in mc.group(1) and "nonsense" not in mc.group(1):
                continue
            out.setdefault(gene.group(1), []).append((f[0], int(f[1]), ref, alt, f[2]))
    return out


def trio_sites(path):
    sites, header, samples = set(), [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                header.append(line)
            elif line.startswith("#CHROM"):
                header.append(line)
                samples = line.rstrip("\n").split("\t")[9:]
            else:
                f = line.split("\t", 3)
                sites.add((f[0], int(f[1])))
    return sites, header, samples


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clinvar", required=True)
    ap.add_argument("--trio", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--truth", required=True)
    a = ap.parse_args()

    cv = clinvar_candidates(a.clinvar)
    sites, header, samples = trio_sites(a.trio)
    order = ["HG002", "HG003", "HG004"]
    assert sorted(samples) == order, samples

    new_recs, truth = [], []
    for name, gene, gts, expect, model in SCENARIOS:
        pool = [c for c in cv.get(gene, []) if (c[0], c[1]) not in sites]
        if len(pool) < len(gts):
            raise SystemExit(f"not enough ClinVar candidates for {gene}")
        chosen = pool[:: max(1, len(pool) // len(gts))][: len(gts)]
        for (chrom, pos, ref, alt, cvid), gt in zip(chosen, gts):
            by_sample = dict(zip(order, gt))
            cols = []
            for s in samples:
                g = by_sample[s]
                ad = {"0/0": "30,0", "0/1": "15,15", "1/1": "0,30"}[g]
                cols.append(f"{g}:{ad}:30:50")
            new_recs.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\tSPIKE=1\tGT:AD:DP:GQ\t" + "\t".join(cols) + "\n")
            truth.append((name, gene, f"{chrom}:{pos}:{ref}:{alt}", expect, model, cvid))
            sites.add((chrom, pos))

    hdr_extra = '##INFO=<ID=SPIKE,Number=0,Type=Flag,Description="Synthetic spiked-in test variant">\n'
    with open(a.out, "w") as out:
        out.writelines(header[:-1])
        out.write(hdr_extra)
        out.write(header[-1])
        with gzip.open(a.trio, "rt") as fh:
            out.writelines(l for l in fh if not l.startswith("#"))
        out.writelines(new_recs)
    with open(a.truth, "w") as t:
        t.write("scenario\tgene\tvariant\texpected\tmodel\tclinvar_id\n")
        for r in truth:
            t.write("\t".join(map(str, r)) + "\n")
    print(f"spiked {len(new_recs)} variants; truth in {a.truth}")


if __name__ == "__main__":
    main()
