#!/usr/bin/env bash
# End-to-end validation on the GIAB Ashkenazi trio (HG002 son, HG003 father, HG004 mother),
# 30x PCR-free NovaSeq WGS, restricted to the primary-immunodeficiency panel genes so it
# runs on a laptop-class machine (4 cores / 16 GB / ~15 GB free disk).
#
# Mirrors the nf-core/raredisease SNV path: DeepVariant -> GLnexus joint calling -> VEP,
# then benchmarks calls against GIAB v4.2.1 truth with hap.py, spikes in known pathogenic
# variants, and runs dnasearch on both the clean and spiked trio.
#
# usage: validation/giab_trio/run.sh [WORKDIR] [STAGE...]
#   stages: slice call joint bench annotate spike report  (default: all, in order)
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
WORK=$(realpath -m "${1:-$REPO/work/giab_trio}")
shift || true
STAGES=${*:-slice call joint bench annotate spike report}
THREADS=${THREADS:-4}
mkdir -p "$WORK"/{ref,bam,dv,annot,bench,reports}

BED="$HERE/pid_panel.grch38.bed"
BAM_BASE=https://storage.googleapis.com/brain-genomics-public/research/sequencing/grch38/bam/novaseq/wgs_pcr_free/30x
REF_URL=https://storage.googleapis.com/deepvariant/case-study-testdata/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio
declare -A TRUTH_DIR=([HG002]=HG002_NA24385_son [HG003]=HG003_NA24149_father [HG004]=HG004_NA24143_mother)
CLINVAR_URL=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
GFF_URL=https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.gff3.gz
GNOMAD=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites

SAMTOOLS=quay.io/biocontainers/samtools:1.21--h96c455f_1
BCFTOOLS=quay.io/biocontainers/bcftools:1.21--h3a4d415_1
DEEPVARIANT=google/deepvariant:1.9.0
GLNEXUS=ghcr.io/dnanexus-rnd/glnexus:v1.4.1
HAPPY=jmcdani20/hap.py:v0.3.12
VEP=ensemblorg/ensembl-vep:release_113.0

# containers reach the network through the same proxy/CA as the host, if one is configured
DOCKER_NET=(--network host)
[[ -n "${HTTPS_PROXY:-}" ]] && DOCKER_NET+=(-e "HTTPS_PROXY=$HTTPS_PROXY" -e "https_proxy=$HTTPS_PROXY")
CA=${CURL_CA_BUNDLE:-${SSL_CERT_FILE:-}}
[[ -n "$CA" && -f "$CA" ]] && DOCKER_NET+=(-v "$CA:/etc/proxy-ca.crt:ro" -e CURL_CA_BUNDLE=/etc/proxy-ca.crt)
run() { local img=$1; shift; docker run --rm -u "$(id -u):$(id -g)" "${DOCKER_NET[@]}" \
          -v "$WORK:$WORK" -v "$HERE:$HERE:ro" -w "$WORK" "$img" "$@"; }
log() { echo "[$(date +%H:%M:%S)] $*"; }

REF="$WORK/ref/GRCh38_no_alt.fna"
stage_ref() {
  [[ -s $REF.fai ]] && return
  log "downloading GRCh38 no-alt analysis set (~3.1 GB)"
  curl -fsSL -o "$REF" "$REF_URL"; curl -fsSL -o "$REF.fai" "$REF_URL.fai"
}

stage_slice() {
  stage_ref
  for s in HG002 HG003 HG004; do
    out="$WORK/bam/$s.pid.bam"
    [[ -s $out.bai ]] && continue
    log "extracting $s reads over panel regions"
    run "$SAMTOOLS" samtools view -@ "$THREADS" -b -M -L "$BED" -o "$out" "$BAM_BASE/$s.novaseq.pcr-free.30x.dedup.grch38.bam"
    run "$SAMTOOLS" samtools index "$out"
  done
  rm -f "$WORK"/*.dedup.grch38.bam.bai  # remote indexes samtools caches in the cwd
}

stage_call() {
  printf 'chrX\t10000\t2781479\nchrX\t155701382\t156030895\nchrY\t10000\t2781479\nchrY\t56887902\t57217415\n' > "$WORK/ref/par.bed"
  for s in HG002 HG003 HG004; do
    out="$WORK/dv/$s.g.vcf.gz"
    [[ -s $out ]] && continue
    hap=()
    [[ $s == HG002 ]] && hap=(--haploid_contigs=chrX,chrY --par_regions_bed="$WORK/ref/par.bed")
    log "DeepVariant $s"
    run "$DEEPVARIANT" /opt/deepvariant/bin/run_deepvariant --model_type=WGS \
      --ref="$REF" --reads="$WORK/bam/$s.pid.bam" --regions="$BED" \
      --output_vcf="$WORK/dv/$s.vcf.gz" --output_gvcf="$out" \
      --intermediate_results_dir="$WORK/dv/tmp_$s" --num_shards="$THREADS" "${hap[@]}"
    rm -rf "$WORK/dv/tmp_$s"
  done
}

stage_joint() {
  [[ -s $WORK/trio.vcf.gz.tbi ]] && return
  log "GLnexus joint genotyping"
  rm -rf "$WORK/GLnexus.DB"
  run "$GLNEXUS" sh -c "glnexus_cli --config DeepVariantWGS --bed $BED --threads $THREADS \
      $WORK/dv/HG002.g.vcf.gz $WORK/dv/HG003.g.vcf.gz $WORK/dv/HG004.g.vcf.gz > $WORK/trio.bcf"
  rm -rf "$WORK/GLnexus.DB"
  run "$BCFTOOLS" sh -c "bcftools norm -m -any -f $REF $WORK/trio.bcf -Oz -o $WORK/trio.vcf.gz && bcftools index -t $WORK/trio.vcf.gz"
}

stage_bench() {
  for s in HG002 HG003 HG004; do
    [[ -s $WORK/bench/$s.summary.csv ]] && continue
    t="$GIAB/${TRUTH_DIR[$s]}/NISTv4.2.1/GRCh38/${s}_GRCh38_1_22_v4.2.1_benchmark"
    log "hap.py $s vs GIAB v4.2.1 (panel regions within the benchmark BED)"
    curl -fsSL -o "$WORK/bench/$s.truth.bed" "${t}_noinconsistent.bed"
    run "$BCFTOOLS" sh -c "bcftools view -R $BED $t.vcf.gz -Oz -o $WORK/bench/$s.truth.vcf.gz && bcftools index -t $WORK/bench/$s.truth.vcf.gz"
    run "$BCFTOOLS" sh -c "bcftools view -s $s -c 1 $WORK/trio.vcf.gz -Oz -o $WORK/bench/$s.query.vcf.gz && bcftools index -t $WORK/bench/$s.query.vcf.gz"
    run "$HAPPY" /opt/hap.py/bin/hap.py "$WORK/bench/$s.truth.vcf.gz" "$WORK/bench/$s.query.vcf.gz" \
      -f "$WORK/bench/$s.truth.bed" -T "$BED" -r "$REF" -o "$WORK/bench/$s" --threads "$THREADS" --engine=vcfeval
  done
  python3 "$HERE/summarize_bench.py" "$WORK/bench" | tee "$WORK/bench/summary.md"
}

stage_annotate() {
  A="$WORK/annot"
  if [[ ! -s $A/clinvar.vcf.gz.tbi ]]; then
    log "ClinVar"
    [[ -s $A/clinvar.raw.vcf.gz ]] || curl -fsSL -o "$A/clinvar.raw.vcf.gz" "$CLINVAR_URL"
    # ClinVar uses bare chromosome names; the reference uses chr-prefixed ones
    for c in $(seq 1 22) X Y MT; do echo "$c chr${c/MT/M}"; done > "$A/chrmap.txt"
    # ClinVar's header has no ##contig lines; bcftools needs an index to know them
    run "$BCFTOOLS" sh -c "tabix -f -p vcf $A/clinvar.raw.vcf.gz && bcftools annotate --rename-chrs $A/chrmap.txt $A/clinvar.raw.vcf.gz -Oz -o $A/clinvar.vcf.gz && bcftools index -t $A/clinvar.vcf.gz"
    rm -f "$A/clinvar.raw.vcf.gz"
  fi
  if [[ ! -s $A/genes.gff3.gz.tbi ]]; then
    log "Ensembl 113 GFF3"
    [[ -s $A/genes.raw.gff3.gz ]] || curl -fsSL -o "$A/genes.raw.gff3.gz" "$GFF_URL"
    gunzip -c "$A/genes.raw.gff3.gz" | grep -v '^#' | awk 'BEGIN{OFS="\t"} $1!~/^(GL|KI)/ {$1=($1=="MT"?"chrM":"chr"$1); print}' \
      | sort -k1,1 -k4,4n -k5,5n -t$'\t' > "$A/genes.gff3"
    run "$SAMTOOLS" sh -c "bgzip -f $A/genes.gff3 && tabix -p gff $A/genes.gff3.gz"
    rm -f "$A/genes.raw.gff3.gz"
  fi
}

# gnomAD v4.1 genome frequencies are fetched remotely by index (-R), only for sites that
# survive the consequence filter, so a run downloads megabytes rather than the ~700 GB release.
fetch_gnomad() {  # fetch_gnomad sites.tsv  (CHROM<TAB>POS); appends to the shared cache
  local sites=$1 A="$WORK/annot"
  touch "$A/gnomad.cache.tsv" "$A/gnomad.queried.tsv"
  LC_ALL=C comm -23 <(LC_ALL=C sort -u "$sites") <(LC_ALL=C sort -u "$A/gnomad.queried.tsv") \
    | sort -k1,1 -k2,2n > "$A/q.new.tsv"
  log "gnomAD v4.1: $(wc -l < "$A/q.new.tsv") new sites to look up"
  for c in $(cut -f1 "$A/q.new.tsv" | sort -u); do
    awk -v c="$c" '$1==c' "$A/q.new.tsv" > "$A/q.$c.tsv"
    run "$BCFTOOLS" bcftools query -R "$A/q.$c.tsv" -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%AF_grpmax\n' \
      "$GNOMAD.$c.vcf.bgz" >> "$A/gnomad.cache.tsv" || { log "WARN gnomAD fetch failed for $c"; continue; }
    cat "$A/q.$c.tsv" >> "$A/gnomad.queried.tsv"
  done
  rm -f "$A"/q.*.tsv
  sort -u "$A/gnomad.cache.tsv" | sort -k1,1V -k2,2n > "$A/gnomad.tsv"
  run "$SAMTOOLS" sh -c "bgzip -f $A/gnomad.tsv && tabix -f -s1 -b2 -e2 $A/gnomad.tsv.gz"
}

annotate_vcf() {  # annotate_vcf in.vcf.gz out.vcf.gz   (ClinVar -> VEP -> gnomAD)
  local in=$1 out=$2 A="$WORK/annot" base=${2%.vcf.gz}
  run "$BCFTOOLS" sh -c "bcftools annotate -a $A/clinvar.vcf.gz -c INFO/CLNSIG,INFO/CLNREVSTAT $in -Oz -o $base.cv.vcf.gz"
  # GFF + FASTA annotation needs no VEP cache (the cache is ~25 GB); --offline would demand one
  run "$VEP" vep --format vcf --vcf --force_overwrite --no_stats --fork "$THREADS" \
    -i "$base.cv.vcf.gz" -o "$base.vep.vcf.gz" --compress_output bgzip \
    --fasta "$REF" --gff "$A/genes.gff3.gz" --symbol --hgvs --canonical --allele_number \
    --fields "Allele,ALLELE_NUM,Consequence,IMPACT,SYMBOL,Gene,Feature,CANONICAL,HGVSc,HGVSp"
  run "$BCFTOOLS" sh -c "bcftools query -i 'INFO/CSQ~\"missense|stop_gained|frameshift|splice_|start_lost|stop_lost|inframe_|protein_altering\" || INFO/CLNSIG!=\".\"' -f '%CHROM\t%POS\n' $base.vep.vcf.gz" > "$base.sites.tsv"
  fetch_gnomad "$base.sites.tsv"
  printf '##INFO=<ID=gnomAD_AF,Number=A,Type=Float,Description="gnomAD v4.1 genomes AF">\n##INFO=<ID=gnomAD_AF_grpmax,Number=A,Type=Float,Description="gnomAD v4.1 genomes grpmax AF">\n' > "$A/gnomad.hdr"
  run "$BCFTOOLS" sh -c "bcftools annotate -a $A/gnomad.tsv.gz -c CHROM,POS,REF,ALT,gnomAD_AF,gnomAD_AF_grpmax -h $A/gnomad.hdr $base.vep.vcf.gz -Oz -o $out"
  rm -f "$base.cv.vcf.gz" "$base.vep.vcf.gz" "$base.sites.tsv"
}

stage_spike() {
  [[ -s $WORK/trio.spiked.vcf.gz ]] && return
  log "spiking known pathogenic variants into the trio"
  python3 -I "$HERE/spike_in.py" --clinvar "$WORK/annot/clinvar.vcf.gz" --trio "$WORK/trio.vcf.gz" \
    --out "$WORK/trio.spiked.unsorted.vcf" --truth "$WORK/spike_truth.tsv"
  run "$BCFTOOLS" sh -c "bcftools sort $WORK/trio.spiked.unsorted.vcf -Oz -o $WORK/trio.spiked.vcf.gz && bcftools index -t $WORK/trio.spiked.vcf.gz"
  rm -f "$WORK/trio.spiked.unsorted.vcf"
}

stage_report() {
  for v in trio trio.spiked; do
    [[ -s $WORK/$v.annotated.vcf.gz ]] || { log "VEP $v"; annotate_vcf "$WORK/$v.vcf.gz" "$WORK/$v.annotated.vcf.gz"; }
  done
  cd "$REPO"
  log "dnasearch: clean trio (expect no Tier 1)"
  python3 -m dnasearch report --case "$HERE/case.toml" --vcf "$WORK/trio.annotated.vcf.gz" --out "$WORK/reports/clean"
  log "dnasearch: spiked trio"
  python3 -m dnasearch report --case "$HERE/case_spiked.toml" --vcf "$WORK/trio.spiked.annotated.vcf.gz" --out "$WORK/reports/spiked"
  python3 "$HERE/check_spike.py" "$WORK/spike_truth.tsv" "$WORK/reports/spiked"
}

for st in $STAGES; do "stage_$st"; done
log "done"
