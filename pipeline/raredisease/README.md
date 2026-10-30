# Running nf-core/raredisease

[nf-core/raredisease](https://nf-co.re/raredisease) is the "call everything" step:

- Alignment (bwa-mem2)
- Small variants with DeepVariant + GLnexus
- Structural variants and copy-number changes (Manta, TIDDIT, CNVnator)
- Repeat expansions (ExpansionHunter + Stranger)
- Mobile elements
- Runs of homozygosity / uniparental disomy
- Mitochondrial variants
- SMN1/SMN2 copy number
- VEP annotation and genmod ranking

`dnasearch` then works on its annotated SNV VCF.

## Hardware

| | singleton 30x | trio |
|---|---|---|
| CPU | 16+ cores | 16-32 cores |
| RAM | 64 GB (32 GB minimum) | 64 GB |
| Disk | ~600 GB working + ~300 GB references | ~1.5 TB working |
| Time | 1-2 days on 16 cores | 2-4 days |

A cloud VM is often easier. On AWS, a `c6i.8xlarge`-class machine with a 2 TB disk works; at
spot prices, expect tens of dollars per genome.

## Steps

1. Install Nextflow (≥ 24.10) and Docker (or Singularity/Apptainer).
2. **Build the references.** This is the main effort:
   - Follow the [raredisease usage docs, "Reference files"](https://nf-co.re/raredisease/3.1.2/docs/usage)
     for each file in `params.yaml`.
   - The Clinical Genomics Stockholm reference bundle used by the pipeline's authors is the
     quickest route.
   - The VEP cache plus CADD alone are roughly 200 GB.
3. Put your FASTQs into `samplesheet.csv`:
   - Copy it from `samplesheet.example.csv`.
   - `sex`: 1=male, 2=female. `phenotype`: 1=unaffected, 2=affected.
   - For a trio, add parents with `paternal_id`/`maternal_id` on the child's rows.
   - You can also start from a vendor BAM using the `bam`/`bai` columns. Only do this if it was
     aligned to the **same** reference FASTA. Otherwise start from FASTQ, or convert BAM to FASTQ
     with `samtools fastq`.
4. Run:

   ```bash
   nextflow run nf-core/raredisease -r 3.1.2 -profile docker -params-file pipeline/raredisease/params.yaml
   ```

5. Triage:

   ```bash
   dnasearch report --case cases/me.toml --vcf /data/results/raredisease/rank_and_filter/<case>_snv_ranked_research.vcf.gz --out results/report
   ```

   Read the frequency and ClinVar fields from the INFO keys your vcfanno config writes. The
   defaults try `GNOMADAF_popmax`, `GNOMADAF` and `CLNSIG`. Override them in the case file's
   `[annotation]` section if yours differ.

Then review the structural-variant, repeat-expansion and mitochondrial outputs
(`*_sv_ranked*.vcf.gz`, `repeat_expansions/`, `*_mito_ranked*.vcf.gz`). Load them into
[Scout](https://github.com/Clinical-Genomics/scout) or
[seqr](https://github.com/broadinstitute/seqr) for interactive filtering.

## Status

The parameter names here were checked against the 3.1.2 schema. This full pipeline has **not**
been run in this repo yet; only the SNV path (DeepVariant → GLnexus → VEP → dnasearch) has been
validated, on GIAB (see `validation/giab_trio`).
