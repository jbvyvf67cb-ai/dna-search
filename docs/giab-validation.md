# Validation on the GIAB Ashkenazi trio

## Setup

- **Data:** HG002 (son), HG003 (father), HG004 (mother) from the Genome in a Bottle project:
  - 30x PCR-free NovaSeq WGS (`brain-genomics-public`), aligned to GRCh38 with bwa-mem;
  - the same kind of data as a consumer 30x genome.
- **Regions:** gene bodies ±5 kb of the 374 green genes in the PanelApp primary immunodeficiency
  panel (v9.112), about 25.8 Mb.
- **Calling:** the SNV path nf-core/raredisease uses:
  - DeepVariant 1.9.0 (WGS model; chrX/Y haploid for HG002 outside the PARs);
  - GLnexus v1.4.1 joint genotyping;
  - `bcftools norm`.
- **Annotation:**
  - ClinVar (GRCh38 weekly VCF);
  - Ensembl VEP 113 with the Ensembl 113 GFF3;
  - gnomAD v4.1 genome allele frequencies (AF, AF_grpmax), queried by index for the
    protein/splice-affecting and ClinVar sites.
- **Phenotype used for ranking:** agammaglobulinaemia, recurrent bacterial infections, recurrent
  respiratory infections. This is hypothetical; the GIAB donors are healthy volunteers. HG002 is
  labelled "affected" only to exercise the trio logic.

Reproduce with `validation/giab_trio/run.sh <workdir>`. It took about 1.5 h on 4 cores and
16 GB RAM, using about 15 GB of disk on top of the container images.

## 1. Variant calling accuracy vs GIAB v4.2.1 truth (hap.py / vcfeval)

Panel regions within each sample's high-confidence benchmark regions; autosomes only (v4.2.1
covers chr1-22).

| sample | type | truth variants | recall | precision | F1 |
|---|---|---|---|---|---|
| HG002 | SNP | 24,598 | 0.9967 | 0.9995 | 0.9981 |
| HG002 | INDEL | 4,205 | 0.9938 | 0.9986 | 0.9962 |
| HG003 | SNP | 24,288 | 0.9954 | 0.9995 | 0.9975 |
| HG003 | INDEL | 4,082 | 0.9927 | 0.9993 | 0.9959 |
| HG004 | SNP | 25,107 | 0.9933 | 0.9993 | 0.9963 |
| HG004 | INDEL | 4,192 | 0.9917 | 0.9976 | 0.9946 |

Trio Mendelian error rate (autosomal, GQ ≥ 20): **0.092%** (38 / 41,344 sites).

The sex check from chrX heterozygosity correctly inferred HG002 as male: 0.00 heterozygous
fraction over 960 non-PAR sites.

**Caveat:** this accuracy holds *inside* GIAB's confident regions. The genes in the report's
"could not assess" list (NCF1, IKBKG, C4A/B, ...) are largely outside those regions, which is
why they are listed.

## 2. Triage on the unmodified trio (healthy family): expect nothing important

| funnel step (HG002's variants) | count |
|---|---|
| all variants in the regions | 35,438 |
| pass quality | 32,738 |
| protein/splice-affecting or ClinVar P/LP | 476 |
| rare (gnomAD) | 13 |
| fit an inheritance model | 12 |
| reported | 2 |

Result: **Tier 1: 0, Tier 2: 0, Tier 3: 2.** Both Tier 3 entries are heterozygous variants
inherited from a healthy parent, flagged "needs reduced penetrance":

- An IFNGR2 missense variant (p.Ile232Val), absent from gnomAD.
- A CHD7 splice-region variant with conflicting ClinVar classifications.

That's the expected outcome for healthy people: a short research list and no diagnostic
candidate. Full report: [`validation/giab_trio/results/giab_trio_clean.report.md`](../validation/giab_trio/results/giab_trio_clean.report.md).

## 3. Spike-in test: known pathogenic variants with designed inheritance

Eight ClinVar "Pathogenic" SNVs (multi-submitter or expert-panel reviewed) were inserted with
chosen genotypes:

| scenario | gene | genotypes (son/father/mother) | expected | result |
|---|---|---|---|---|
| X-linked (XLA) | BTK | 1/1, 0/0, 0/1 | Tier 1 x-linked hemizygous | **PASS** |
| de novo dominant (STAT3-HIES) | STAT3 | 0/1, 0/0, 0/0 | Tier 1 de novo | **PASS** |
| compound het, one from each parent (RAG1 SCID) | RAG1 ×2 | 0/1,0/1,0/0 + 0/1,0/0,0/1 | Tier 1 compound het | **PASS** |
| **negative control:** two variants both from the father (cis) | RAG2 ×2 | 0/1,0/1,0/0 ×2 | not Tier 1/2 (carrier only) | **PASS** |
| single recessive allele | DOCK8 | 0/1, 0/0, 0/1 | carrier list only | **PASS** |
| incidental finding | BRCA2 | 0/1, 0/1, 0/0 | ACMG SF section (opt-in) | **PASS** |

Report: [`validation/giab_trio/results/giab_trio_spiked.report.md`](../validation/giab_trio/results/giab_trio_spiked.report.md).

## Bugs the validation caught (all fixed, with regression tests)

1. **Neighbouring genes borrowed evidence.** RAG1 and RAG2 are adjacent, so a RAG1 variant's
   "downstream of RAG2" annotation produced a fake RAG2 compound het, and the BTK variant also
   counted for TIMM8A. Now only in-gene consequences count.
2. **ClinVar evidence attached to the wrong gene.** A RAG2 ClinVar variant also lies in an RAG1
   intron, and its "Pathogenic" label was counted for RAG1. ClinVar evidence now only counts
   for the gene ClinVar names (`GENEINFO`).
3. **Healthy parents weren't used as evidence against.** The first clean-trio run had a Tier 1
   PARN frameshift for which the healthy father is **homozygous**, which rules out a dominant
   disease. Now:
   - an unaffected homozygous parent removes dominant models;
   - an unaffected parent who is homozygous for, or carries both, variants of a pair removes
     the compound-het model.
4. **Population frequencies were silently missing.** A bcftools filter quirk (no `|`
   alternation in `~`) meant most coding variants never got a gnomAD lookup. They were treated
   as rare, which inflated Tier 2/3. `dnasearch` now reports the number of protein-affecting
   variants without a frequency, and warns when it exceeds 30%.

## What this does not validate

- Structural variants, repeat expansions, mitochondrial variants, mosaicism, and
  paralog-confounded genes. These come from nf-core/raredisease and the extras, which were not
  run here.
- Running the full nf-core/raredisease pipeline. Only its SNV path was reproduced. Its INFO
  field names for gnomAD/ClinVar may need mapping in the case file's `[annotation]` section.
- Real-patient sensitivity. Spike-ins test the logic, not the read-level detectability of each
  variant.
- SpliceAI/CADD scores. They weren't annotated here, so splice-only variants rely on VEP
  consequence alone.
