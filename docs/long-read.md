# Getting a long-read genome (US, as of Oct 2026)

**Summary: no current direct-to-consumer long-read genome was found.** The realistic routes are
a physician-ordered clinical test, a research or core-facility service, or a research study.
Prices change, and several of the facts below come from a single source. Confirm with the
provider before ordering.

## Why long reads here

Long reads resolve what 30x short reads can't:

- Paralog/pseudogene genes (NCF1, IKBKG, SMN1, CYP21A2, C4A/B)
- Structural variants and repeat expansions
- Phasing: which parental copy carries each variant
- Methylation, from the same run

## 1. Direct-to-consumer

| Provider | Status |
|---|---|
| Dante Labs | Current site lists only short-read 30x WGS ([dantelabs.com](https://dantelabs.com/whole-genome-sequencing)). The 2019 ONT "GenomeL" (30x, $999) and the hybrid "Whole GenomeH" appear discontinued ([PR Newswire 2019](https://www.prnewswire.com/news-releases/dante-labs-launches-genomel-the-first-commercial-long-reads-human-whole-genome-sequencing-300826213.html)). Customer-reported turnaround is very long. |
| Nucleus, Sequencing.com, MyHeritage, Human Longevity | No long-read product found |
| Full Genomes | A review site mentions "long-read options"; not verified on the vendor's own site |

## 2. Clinical (physician-ordered)

| Lab | What | Notes |
|---|---|---|
| **Azenta / GENEWIZ** | Clinical long-read WGS on PacBio Revio HiFi; CLIA/CAP; launched Oct 2024 ([Azenta](https://investors.azenta.com/2024-10-30-Azenta-Obtains-Regulatory-Approval-for-Clinical-Long-Read-Whole-Genome-Sequencing-Test)) | Most direct clinical route to a full HiFi genome. Price, ordering and raw-data policy not published; request a quote |
| **GeneDx** | Early-access long-read WGS "by year-end" 2026 for eligible patients selected by phenotype ([announcement](https://www.stocktitan.net/news/WGS/gene-dx-expands-diagnostic-leadership-with-long-read-genome-bwjg58pd7oru.html)) | Platform and price not disclosed |
| **Baylor Genetics** | Since Mar 2026, long reads and optical mapping are reflexed automatically within clinical WGS when a qualifying variant needs them ([CLP](https://clpmag.com/diagnostic-technologies/molecular-diagnostics/baylor-genetics-adds-ogm-lrs-whole-genome-sequencing/)) | Can't be ordered as a standalone test |
| **Children's Mercy Kansas City** | HiFi "5-base" genome (sequence + methylation), in clinical use since 2023 | Pediatric/inpatient focus; outside or adult orders unconfirmed |

For any clinical lab, ask **before ordering** whether you can get the raw FASTQ/BAM. A HIPAA
records request usually works, but labs differ.

## 3. Research / service labs (you pay; results are research-grade)

| Provider | Approx. price | Notes |
|---|---|---|
| Broad Clinical Labs | from $599 (PacBio Revio, 10-40x) ([site](https://broadclinicallabs.org/whole-genome-sequencing/)) | Delivers BAM/FASTQ; every project goes through onboarding/regulatory review; acceptance of individuals unclear |
| Novogene | ~<$1,000 for 20x HiFi (2025 promotion) | Individual acceptance unclear |
| UCSD IGM core | ~$3,500 for ~30x (extraction + library + Revio cell, Feb 2025 external rates) ([rates](https://gcore.ucsd.edu/long-read-prices2)) | Usually needs an institutional account |
| Johns Hopkins core | ~$2,000 for 30x (undated) ([site](https://igc.jhmi.edu/service/long-read-sequencing/)) | Usually needs an institutional account |

Cores generally require a PI or institutional account. Going through a physician-researcher is
often the practical way in.

## 4. Research programs

- **NIH Undiagnosed Diseases Network**:
  - You apply through the UDN Gateway with a physician referral letter.
  - It is aimed at people evaluated by several specialists with no diagnosis.
  - Long reads are used at some sites, but not routinely.
- Disease-specific long-read studies come and go. Search ClinicalTrials.gov for "long-read sequencing" plus your phenotype.

## Practical caveats

- **Blood beats saliva**:
  - Saliva HiFi reads are shorter (~15 kb vs ~17 kb mean).
  - Saliva DNA yields are lower, and it contains microbial DNA.
  - Request EDTA blood with high-molecular-weight extraction.
- **Coverage**: ~30x HiFi is the clinical-study norm. 20x is workable for structural variants
  but weaker for heterozygous small variants.
- **Ask for**:
  - Aligned or unaligned BAM with methylation tags (MM/ML).
  - ONT POD5 signal files, if Nanopore, so you can re-basecall later.
- **Analysis**: use the PacBio WGS pipeline (pbmm2, DeepVariant, pbsv/sawfish, TRGT,
  **Paraphase** for NCF1/SMN1/etc.) or the ONT `wf-human-variation`. Standard short-read
  pipelines will not resolve paralog genes.

Research compiled Oct 2026 from vendor pages and news sources linked above. Treat prices as
indicative.
