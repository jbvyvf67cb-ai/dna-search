# dna-search

A personal-genome workup for an undiagnosed condition, focused on primary immunodeficiency
(inborn errors of immunity). It works in two steps:

1. **Detect everything** in a 30x whole genome. That means SNVs/indels, structural variants,
   repeat expansions, mitochondrial variants, ROH/UPD and more, using
   [nf-core/raredisease](https://nf-co.re/raredisease) plus the add-ons in `pipeline/extras`.
2. **Narrow it down** with `dnasearch`:
   - quality, consequence and population-frequency filters;
   - trio-aware inheritance models;
   - the Genomics England PanelApp primary-immunodeficiency panel;
   - genome-wide Mendelian disease genes from HPO;
   - HPO phenotype ranking;
   - an explicit list of what the data **could not** assess.

> **Research use only.** This produces a triage list for expert review, not a diagnosis. Any
> finding must be confirmed in a CLIA-certified lab and interpreted by a clinical geneticist or
> immunologist. Basic immune labs (immunoglobulins, vaccine titres, lymphocyte subsets,
> complement) should come first; most CVID has no single-gene cause.

## Layout

```
dnasearch/              triage + report package (Python ≥3.11, stdlib only)
cases/example.toml      case file template: pedigree, HPO terms, panels, options
resources/
  panels/               PanelApp primary immunodeficiency panel snapshot (refresh: fetch-resources)
  acmg_sf.tsv           ACMG SF v3.3 secondary-findings genes and reporting rules (opt-in)
  difficult_genes.tsv   genes 30x short reads resolve poorly (NCF1, IKBKG, C4A/B, 22q11.2, ...)
  hpo/                  HPO ontology + annotations (downloaded, git-ignored)
pipeline/raredisease/   params + samplesheet for the "detect everything" pipeline
pipeline/extras/        HLA/KIR, telomere length, somatic, non-human reads, PGx, ...
validation/giab_trio/   end-to-end test on the GIAB Ashkenazi trio (results in docs/giab-validation.md)
docs/                   survey design, long-read options, validation results
```

## Quick start

```bash
pip install -e '.[test]'
dnasearch fetch-resources            # HPO files + refresh the PanelApp panel
cp cases/example.toml cases/me.toml  # fill in sample id, sex, parents, HPO terms
dnasearch report --case cases/me.toml --vcf <annotated.vcf.gz> --out results/me
```

The output is `results/me/<case>.report.md` (for reading) plus `.findings.tsv` and `.json`.

`--vcf` must be a VEP-annotated VCF, such as raredisease's ranked SNV VCF. Allele frequency,
ClinVar and SpliceAI values are read from INFO or CSQ fields; the field lists can be changed in
the case file's `[annotation]` section.

## Validation

`validation/giab_trio/run.sh` does the following with Docker:

- pulls the GIAB HG002/HG003/HG004 30x reads over the 374 panel genes;
- calls with DeepVariant + GLnexus (the same SNV path raredisease uses);
- benchmarks against GIAB v4.2.1 truth with hap.py;
- annotates with VEP + ClinVar + gnomAD v4.1;
- spikes in known pathogenic variants (X-linked, de novo, compound het, a cis negative control,
  carrier, secondary finding);
- checks that `dnasearch` reports each one the way it should.

See [docs/giab-validation.md](docs/giab-validation.md) for results.

```bash
python3 -m pytest          # unit tests
validation/giab_trio/run.sh work/giab_trio   # ~1-2 h on 4 cores; ~15 GB disk
```

## Further reading

- [docs/survey.md](docs/survey.md): what a comprehensive survey covers and how the funnel narrows it
- [docs/long-read.md](docs/long-read.md): how to get a long-read genome in the US
- [pipeline/raredisease/README.md](pipeline/raredisease/README.md): running the full pipeline on your data
