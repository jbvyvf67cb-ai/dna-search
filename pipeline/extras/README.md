# Add-on analyses outside nf-core/raredisease

These run on the aligned CRAM/BAM from raredisease. They are listed with the tool and the
question each one answers. **None has been validated in this repo yet.**

| Analysis | Why (immunology angle) | Tool | Notes |
|---|---|---|---|
| HLA class I/II typing | autoimmunity, drug hypersensitivity, transplant | [T1K](https://github.com/mourisl/T1K) or HLA-HD | 30x WGS gives 2-field types reliably for class I |
| KIR genotyping | NK-cell function | T1K (KIR mode) | copy number less reliable than gene presence |
| Telomere length | telomere biology disorders (dyskeratosis congenita) cause immunodeficiency and marrow failure | [TelSeq](https://github.com/zd1/telseq) | relative measure; compare with age-matched data or clinical flow-FISH |
| Non-human reads | inherited chromosomally integrated HHV-6, EBV/CMV load, contamination | `samtools view -f 12` + Kraken2 | ciHHV-6 gives ~1 viral genome copy per cell |
| Somatic / mosaic | VEXAS (UBA1), clonal haematopoiesis, somatic autoinflammatory (NLRP3, TLR8) | GATK Mutect2 tumour-only mode with a panel of normals | 30x detects clones above ~10-15% VAF only |
| Paralog genes | NCF1 (CGD), IKBKG, CYP21A2, PMS2, GBA1 | [Parascopy](https://github.com/tprodanov/parascopy) | short-read; long reads + Paraphase are better |
| Pharmacogenomics | drug choice and dosing (e.g. immunosuppressants: TPMT/NUDT15 for azathioprine) | [PharmCAT](https://pharmcat.org) + Cyrius (CYP2D6) | |
| Polygenic scores | context only | [pgsc_calc](https://github.com/PGScatalog/pgsc_calc) | ancestry-matched comparison required |
