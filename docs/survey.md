# Comprehensive genome survey for an undiagnosed case

The approach used by undiagnosed-disease programs:

1. Detect **every** kind of variation the data allows. That produces millions of mostly
   meaningless results.
2. Narrow them in stages, using population frequency, inheritance and the patient's phenotype.

This repo implements the narrowing for small variants (`dnasearch`). It delegates detection to
nf-core/raredisease and the add-ons in `pipeline/extras`.

## Layer 1: what a 30x short-read genome can show

| Variant class | Why it matters | Where it's handled |
|---|---|---|
| Coding SNVs and small indels | Most diagnoses | raredisease → **dnasearch** |
| Splice-altering variants, including deep intronic | Missed by panels | VEP + SpliceAI (raredisease) → dnasearch `spliceai_min` |
| Non-coding regulatory variants | Mostly noise; a growing share of diagnoses | Genomiser (not yet integrated) |
| CNVs / structural variants | Whole-gene deletions and duplications | raredisease (Manta, TIDDIT, CNVnator) |
| Repeat expansions | Neurological and multisystem disease | raredisease (ExpansionHunter + Stranger) |
| Mitochondrial DNA, including heteroplasmy | Multisystem disease | raredisease mito workflow |
| Somatic / mosaic in blood | VEXAS (UBA1), clonal haematopoiesis | extras: Mutect2 (30x detects clones above ~10-15% of cells only) |
| Runs of homozygosity / UPD | Hidden recessive disease | raredisease (UPD + ROH) |
| Paralog genes (NCF1, IKBKG, SMN1, CYP21A2) | Standard callers fail | SMNCopyNumberCaller (raredisease); Parascopy (extras); long reads |
| HLA / KIR types | Autoimmunity, drug hypersensitivity | extras: T1K |
| Telomere length | Telomere disorders → immunodeficiency | extras: TelSeq |
| Sex-chromosome aneuploidy / mosaic loss of Y | Often undiagnosed | Coverage QC (raredisease, ngs-bits) |
| Non-human reads | Inherited integrated HHV-6, viral load | extras: Kraken2 |
| Pharmacogenomics, polygenic scores | Context, not diagnosis | extras: PharmCAT, pgsc_calc |

## Layer 2: data beyond one genome (rough order of added value)

1. **Phenotype as HPO terms.** The most important input to the narrowing, and free.
2. **Parents' genomes (trio).** They reveal de novo variants and resolve phase.
3. **RNA-seq** (blood or fibroblasts). Detects aberrant splicing or expression; adds roughly
   10-15% to diagnostic yield. Pipeline: DROP.
4. **Long reads** (PacBio HiFi / ONT). Resolve paralog genes, structural variants, repeats,
   phase and methylation. See [long-read.md](long-read.md).
5. **Methylation episignatures** (clinical EpiSign).
6. **Matchmaking.** Matchmaker Exchange and GeneMatcher, through a clinician or researcher.
7. **Periodic reanalysis.** Roughly 250 new gene-disease associations are published each year.

## The funnel (dnasearch)

```
all proband variants
  → pass quality (FILTER, GQ, depth)
  → affect a protein or splicing, or ClinVar P/LP
  → rare (gnomAD; stricter cut-off for dominant models than recessive)
  → fit an inheritance model (de novo, dominant, homozygous, compound het, X-linked)
  → tier:
      1  panel gene, mode of inheritance fits, strong evidence (ClinVar P/LP, LoF, de novo)
      2  panel gene with weaker evidence, or other Mendelian disease gene with strong evidence
      3  candidate: gene with no known disease (strong, informative model), or weak evidence
         with a good phenotype match
  + opt-in: ACMG SF v3.3 secondary findings; carrier status
  + always: list of what could NOT be assessed
```

Within each tier, results are sorted by HPO phenotype similarity: Resnik best-match average
against HPO gene annotations.

## Pitfalls

- **A broad search always finds plausible-looking candidates.** Each person carries roughly 100
  loss-of-function variants and dozens of ClinVar "pathogenic" calls that are wrong or
  harmless in one copy. Require that the inheritance fits, the gene-disease link is established,
  and the phenotype matches.
- **Incidental findings.** Decide on `secondary_findings` before looking at results.
- **Yield.** Genome sequencing solves roughly 25-40% of undiagnosed cases. A negative result
  does not rule out a genetic cause.
- **Confirmation.** Every finding needs confirmation by a clinical lab and expert review.
