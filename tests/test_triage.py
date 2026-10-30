import textwrap

import pytest

from dnasearch import inheritance as inh
from dnasearch import report, triage
from dnasearch.vcf import HET, HOM_ALT, HOM_REF, Call, VcfReader, vep_allele

CSQ_FMT = "Allele|ALLELE_NUM|Consequence|IMPACT|SYMBOL|Gene|Feature|CANONICAL|HGVSc|HGVSp"


def csq(gene, cons, allele="T", num="1"):
    return f"{allele}|{num}|{cons}|X|{gene}|ENSG|ENST|YES|c.1A>T|p.X1Y"


def rec(chrom, pos, ref, alt, gts, info="", flt="PASS"):
    samples = "\t".join(f"{g}:20,20:40:50" for g in gts)
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t50\t{flt}\t{info or '.'}\tGT:AD:DP:GQ\t{samples}"


HEADER = textwrap.dedent(f"""\
    ##fileformat=VCFv4.2
    ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {CSQ_FMT}">
    ##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="af">
    ##INFO=<ID=CLNSIG,Number=.,Type=String,Description="clinvar">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tKID\tDAD\tMOM
    """)

#                 KID    DAD    MOM
RECORDS = [
    # AD gene, de novo missense -> Tier 1
    rec("chr1", 100, "A", "T", ["0/1", "0/0", "0/0"], f"CSQ={csq('ADGENE', 'missense_variant')}"),
    # AR gene, trans compound het (one from each parent), LoF + ClinVar P -> Tier 1
    rec("chr2", 100, "A", "T", ["0/1", "0/1", "0/0"], f"CSQ={csq('ARGENE', 'stop_gained')}"),
    rec("chr2", 200, "A", "T", ["0/1", "0/0", "0/1"], f"CLNSIG=Pathogenic;CSQ={csq('ARGENE', 'missense_variant')}"),
    # AR gene, both from father (cis) -> not a compound het
    rec("chr3", 100, "A", "T", ["0/1", "0/1", "0/0"], f"CSQ={csq('ARGENE2', 'stop_gained')}"),
    rec("chr3", 200, "A", "T", ["0/1", "0/1", "0/0"], f"CSQ={csq('ARGENE2', 'frameshift_variant')}"),
    # X-linked gene, hemizygous son, carrier mother -> Tier 1
    rec("chrX", 50_000_000, "A", "T", ["1/1", "0/0", "0/1"], f"CSQ={csq('XGENE', 'stop_gained')}"),
    # AR homozygous LoF, carrier parents -> Tier 1
    rec("chr4", 100, "A", "T", ["1/1", "0/1", "0/1"], f"CSQ={csq('ARGENE3', 'frameshift_variant')}"),
    # common de novo-looking variant: filtered by frequency
    rec("chr5", 100, "A", "T", ["0/1", "0/0", "0/0"], f"gnomAD_AF=0.2;CSQ={csq('ADGENE2', 'stop_gained')}"),
    # ClinVar benign LoF: dropped
    rec("chr6", 100, "A", "T", ["0/1", "0/0", "0/0"], f"CLNSIG=Benign;CSQ={csq('ADGENE3', 'stop_gained')}"),
    # synonymous: dropped at effect step
    rec("chr7", 100, "A", "T", ["0/1", "0/0", "0/0"], f"CSQ={csq('ADGENE', 'synonymous_variant')}"),
    # failing FILTER
    rec("chr8", 100, "A", "T", ["0/1", "0/0", "0/0"], f"CSQ={csq('ADGENE', 'stop_gained')}", flt="LowQual"),
    # AD panel gene het inherited from unaffected father, LoF -> Tier 2 (downgraded, reduced penetrance)
    rec("chr9", 100, "A", "T", ["0/1", "0/1", "0/0"], f"CSQ={csq('ADGENE4', 'stop_gained')}"),
    # Mendelian error (kid hom-alt, both parents hom-ref)
    rec("chr10", 100, "A", "T", ["1/1", "0/0", "0/0"], f"CSQ={csq('NOVEL', 'intron_variant')}"),
]

PANEL = textwrap.dedent("""\
    symbol\tconfidence\tmode_of_inheritance
    ADGENE\t3\tMONOALLELIC, autosomal or pseudoautosomal, NOT imprinted
    ADGENE2\t3\tMONOALLELIC, autosomal or pseudoautosomal, NOT imprinted
    ADGENE3\t3\tMONOALLELIC, autosomal or pseudoautosomal, NOT imprinted
    ADGENE4\t3\tMONOALLELIC, autosomal or pseudoautosomal, NOT imprinted
    ARGENE\t3\tBIALLELIC, autosomal or pseudoautosomal
    ARGENE2\t3\tBIALLELIC, autosomal or pseudoautosomal
    ARGENE3\t3\tBIALLELIC, autosomal or pseudoautosomal
    XGENE\t3\tX-LINKED: hemizygous mutation in males, biallelic mutations in females
    """)

CASE = textwrap.dedent("""\
    [case]
    id = "t"
    proband = "KID"
    [[samples]]
    id = "KID"
    sex = "male"
    affected = true
    father = "DAD"
    mother = "MOM"
    [[samples]]
    id = "DAD"
    sex = "male"
    affected = false
    [[samples]]
    id = "MOM"
    sex = "female"
    affected = false
    [resources]
    panels = [{ name = "test", path = "panel.tsv" }]
    hpo_dir = "nonexistent"
    [report]
    carrier_status = true
    """)


@pytest.fixture
def result(tmp_path):
    (tmp_path / "t.vcf").write_text(HEADER + "\n".join(RECORDS) + "\n")
    (tmp_path / "panel.tsv").write_text(PANEL)
    (tmp_path / "case.toml").write_text(CASE)
    case = triage.load_case(str(tmp_path / "case.toml"))
    return case, triage.run(case, str(tmp_path / "t.vcf")), tmp_path


def by_gene(res):
    out = {}
    for f in res.findings:
        out.setdefault(f.gene, []).append((f.tier, f.model))
    return out


def test_tiers(result):
    _, res, _ = result
    g = by_gene(res)
    assert ("1", inh.DE_NOVO) in g["ADGENE"]
    assert ("1", inh.COMPHET) in g["ARGENE"]
    assert ("1", inh.XL_HEMI) in g["XGENE"]
    assert ("1", inh.AR_HOM) in g["ARGENE3"]
    assert g["ADGENE4"] == [("2", inh.AD_INHERITED)]
    assert "ARGENE2" not in g  # cis pair is not a compound het
    assert "ADGENE2" not in g  # common
    assert "ADGENE3" not in g  # ClinVar benign


def test_funnel_and_qc(result):
    _, res, _ = result
    funnel = dict(res.funnel)
    assert funnel["proband variants"] == 13
    assert funnel["pass quality"] == 12
    assert res.qc["mendelian_errors"] == 5  # 4 PASS de novo-pattern sites + 1 hom-alt error


def test_carriers_exclude_solved_genes(result):
    _, res, _ = result
    carriers = {f.gene for f in res.carriers}
    assert "ARGENE2" in carriers
    assert "ARGENE" not in carriers  # explained as compound het


def test_report_writes(result):
    case, res, tmp = result
    paths = report.write(case, res, str(tmp / "out"))
    md = open(paths["markdown"]).read()
    assert "Tier 1" in md and "ADGENE" in md
    assert "Not evaluated" in md  # secondary findings disabled


def test_vep_allele():
    assert vep_allele("A", "T", ["T"]) == "T"
    assert vep_allele("AT", "A", ["A"]) == "-"
    assert vep_allele("A", "AT", ["AT"]) == "T"
    assert vep_allele("A", "AT", ["AT", "G"]) == "AT"  # mixed: VEP keeps full alleles


def test_multiallelic_split(tmp_path):
    p = tmp_path / "m.vcf"
    p.write_text(HEADER + rec("chr1", 5, "A", "T,G", ["1/2", "0/1", "0/2"]) + "\n")
    vs = list(VcfReader(str(p)))
    assert [v.alt for v in vs] == ["T", "G"]
    assert vs[0].calls["KID"].gt == HET and vs[0].calls["MOM"].gt == HOM_REF
    assert vs[1].calls["MOM"].gt == HET


def test_mendelian_error():
    c = lambda g: Call(gt=g)  # noqa: E731
    assert inh.mendelian_error(c(HOM_ALT), c(HOM_REF), c(HET)) is True
    assert inh.mendelian_error(c(HET), c(HOM_REF), c(HET)) is False
    assert inh.mendelian_error(c(None), c(HOM_REF), c(HET)) is None


def test_secondary_findings_rules(tmp_path):
    recs = [
        # single P/LP in an AD gene -> reported
        rec("chr13", 100, "A", "T", ["0/1", "0/1", "0/0"], f"CLNSIG=Pathogenic;CSQ={csq('BRCA2', 'stop_gained')}"),
        # single P/LP in a 2-allele gene -> not reported
        rec("chr1", 300, "A", "T", ["0/1", "0/1", "0/0"], f"CLNSIG=Pathogenic;CSQ={csq('MUTYH', 'missense_variant')}"),
        # TTN missense P -> not reported (truncating only)
        rec("chr2", 500, "A", "T", ["0/1", "0/1", "0/0"], f"CLNSIG=Pathogenic;CSQ={csq('TTN', 'missense_variant')}"),
    ]
    (tmp_path / "t.vcf").write_text(HEADER + "\n".join(recs) + "\n")
    (tmp_path / "panel.tsv").write_text(PANEL)
    (tmp_path / "sf.tsv").write_text(
        "symbol\tcategory\tinheritance\treport_rule\tmin_alleles\tspecial\n"
        "BRCA2\tCancer\tAD\tAll P and LP\t1\t\n"
        "MUTYH\tCancer\tAR\tP and LP (2 variants)\t2\t\n"
        "TTN\tCardio\tAD\ttruncating only\t1\tlof_only\n")
    (tmp_path / "case.toml").write_text(
        CASE.replace('hpo_dir = "nonexistent"', 'hpo_dir = "nonexistent"\nacmg_sf = "sf.tsv"')
        .replace("[report]", "[report]\nsecondary_findings = true"))
    case = triage.load_case(str(tmp_path / "case.toml"))
    res = triage.run(case, str(tmp_path / "t.vcf"))
    assert [f.gene for f in res.secondary] == ["BRCA2"]
