#!/bin/bash
set -euo pipefail

# Samples and control 
samples=(
  "SRR21423920"
)
default_input="SRR21423921"

# Project structure
PROJECT_DIR="ppath_chip"
BAM_DIR="${PROJECT_DIR}/out/bams"
MACS2_OUT="${PROJECT_DIR}/out/macs2_output"
HOMER_OUT="${PROJECT_DIR}/out/homer_output"
MEME_OUT="${PROJECT_DIR}/out/meme_output"

QUALITY_OUT="${PROJECT_DIR}/out/qc"
FILTERED_BAM_DIR="${PROJECT_DIR}/out/filtered_bams"
TRACKS_OUT="${PROJECT_DIR}/out/tracks"
FRIP_OUT="${PROJECT_DIR}/out/frip_output"
IDR_OUT="${PROJECT_DIR}/out/idr_output"

# Reference files
REF_GENOME_FA="${PROJECT_DIR}/in/hg38/genome.fa"
REF_GTF="${PROJECT_DIR}/in/hg38/hg38.refGene.gtf"
BLACKLIST_BED="${PROJECT_DIR}/in/hg38/hg38-blacklist.v2.bed"
REF_chr_sizes="${PROJECT_DIR}/in/hg38/chrom.sizes"

cd "${PROJECT_DIR}"

for sample in "${samples[@]}"; do
    echo "[ $(date '+%Y-%m-%d %H:%M:%S') ] Processing sample: ${sample}"

    # raw BAMs
    bam_t="${BAM_DIR}/${sample}/${sample}.bam"
    bam_c="${BAM_DIR}/${default_input}/${default_input}.bam"

    # Check inputs
    if [[ ! -f "${bam_t}" ]]; then
        echo "[ WARNING ] Treatment BAM not found for '${sample}': ${bam_t}. Skipping."
        continue
    fi
    if [[ ! -f "${bam_c}" ]]; then
        echo "[ ERROR ] Control BAM not found: ${bam_c}. Exiting."
        exit 1
    fi

    # prepare outputs
    mkdir -p \
      "${MACS2_OUT}/${sample}" \
      "${HOMER_OUT}/${sample}" \
      "${MEME_OUT}/${sample}" \
      "${QUALITY_OUT}/${sample}" \
      "${FILTERED_BAM_DIR}/${sample}" \
      "${TRACKS_OUT}/${sample}" \
      "${FRIP_OUT}/${sample}" \
      "${IDR_OUT}/${sample}"

    #
    # 0) Pre-peak QC: flagstat, idxstats, insert-size, cross-correlation
    #

    samtools flagstat "${bam_t}" > "${QUALITY_OUT}/${sample}/${sample}.flagstat"
    samtools idxstats "${bam_t}" > "${QUALITY_OUT}/${sample}/${sample}.idxstats"

    module load picard/2.21.6-Java-11.0.2
    java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics \
        INPUT="${bam_t}" \
        OUTPUT="${QUALITY_OUT}/${sample}/${sample}_insert_size_metrics.txt" \
        HISTOGRAM_FILE="${QUALITY_OUT}/${sample}/${sample}_insert_size_hist.pdf" \
        VALIDATION_STRINGENCY=LENIENT

    # 1) Remove chrM, mark duplicates, filter MAPQ & proper pairs

    # remove mitochondrial reads
    samtools view -h "${bam_t}" | grep -v "^.*chrM.*$" \
      | samtools sort -O BAM -o "${FILTERED_BAM_DIR}/${sample}/${sample}_noChrM.bam"
    samtools index "${FILTERED_BAM_DIR}/${sample}/${sample}_noChrM.bam"

    module load picard/2.21.6-Java-11.0.2
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        INPUT="${FILTERED_BAM_DIR}/${sample}/${sample}_noChrM.bam" \
        OUTPUT="${FILTERED_BAM_DIR}/${sample}/${sample}_dedup.bam" \
        METRICS_FILE="${QUALITY_OUT}/${sample}/${sample}_dup_metrics.txt" \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT

    # keep only proper pairs MAPQ>=30
    samtools view -@ 10 -b -f 2 -q 30 \
        "${FILTERED_BAM_DIR}/${sample}/${sample}_dedup.bam" \
      | samtools sort -O BAM \
            -o "${FILTERED_BAM_DIR}/${sample}/${sample}_filtered.bam"
    samtools index "${FILTERED_BAM_DIR}/${sample}/${sample}_filtered.bam"

    # repeat filtering for control BAM
    samtools view -h "${bam_c}" | grep -v "^.*chrM.*$" \
      | samtools sort -O BAM -o "${FILTERED_BAM_DIR}/${sample}/${default_input}_noChrM.bam"
    samtools index "${FILTERED_BAM_DIR}/${sample}/${default_input}_noChrM.bam"
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        INPUT="${FILTERED_BAM_DIR}/${sample}/${default_input}_noChrM.bam" \
        OUTPUT="${FILTERED_BAM_DIR}/${sample}/${default_input}_dedup.bam" \
        METRICS_FILE="${QUALITY_OUT}/${sample}/${default_input}_dup_metrics.txt" \
        REMOVE_DUPLICATES=true \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT
    samtools view -@ 10 -b -f 2 -q 30 \
        "${FILTERED_BAM_DIR}/${sample}/${default_input}_dedup.bam" \
      | samtools sort -O BAM \
            -o "${FILTERED_BAM_DIR}/${sample}/${default_input}_filtered.bam"
    samtools index "${FILTERED_BAM_DIR}/${sample}/${default_input}_filtered.bam"

    # point to filtered BAMs for peak calling
    bam_t="${FILTERED_BAM_DIR}/${sample}/${sample}_filtered.bam"
    bam_c="${FILTERED_BAM_DIR}/${sample}/${default_input}_filtered.bam"


    # 2) MACS2 peak calling
    macs2 callpeak \
        -t "${bam_t}" \
        -c "${bam_c}" \
        -f BAM \
        --call-summits \
        -g hs \
        -B --SPMR \
        --outdir "${MACS2_OUT}/${sample}" \
        -n "${sample}" \
        --cutoff-analysis

    # 3) Generate browser tracks (bedGraph → bigWig)
    sort -k1,1 -k2,2n \
      "${MACS2_OUT}/${sample}/${sample}_treat_pileup.bdg" \
      > "${TRACKS_OUT}/${sample}/${sample}_pileup.sorted.bdg"
    ./src/bedGraphToBigWig \
      "${TRACKS_OUT}/${sample}/${sample}_pileup.sorted.bdg" \
      "${REF_chr_sizes}" \
      "${TRACKS_OUT}/${sample}/${sample}_pileup.bw"

    # 4) Filter peaks against blacklist
    bedtools intersect -v \
        -a "${MACS2_OUT}/${sample}/${sample}_peaks.narrowPeak" \
        -b "${BLACKLIST_BED}" \
      > "${MACS2_OUT}/${sample}/${sample}_noblacklist.narrowPeak"

    # 5) FRiP score
    # create SAF for peaks
    awk 'BEGIN{FS=OFS="\t"} \
        { if (NF==10) print $4, $1, $2+1, $3, "+"; }' \
        "${MACS2_OUT}/${sample}/${sample}_noblacklist.narrowPeak" \
      > "${FRIP_OUT}/${sample}/${sample}.saf"
    featureCounts \
        -p \
        -F SAF \
        -a "${FRIP_OUT}/${sample}/${sample}.saf" \
        -o "${FRIP_OUT}/${sample}/${sample}_readCountInPeaks.txt" \
        "${bam_t}"

    # 6) HOMER annotation & motif finding
    cut -f1-3 "${MACS2_OUT}/${sample}/${sample}_noblacklist.narrowPeak" \
        > "${HOMER_OUT}/${sample}/${sample}.bed"
    annotatePeaks.pl \
        "${HOMER_OUT}/${sample}/${sample}.bed" \
        hg38 \
        > "${HOMER_OUT}/${sample}/${sample}_annotatedPeaks.txt"
    findMotifsGenome.pl \
        "${HOMER_OUT}/${sample}/${sample}.bed" \
        hg38 \
        "${HOMER_OUT}/${sample}" \
        -preparsedDir "${HOMER_OUT}/${sample}/preparsed"

    # 7) MEME-ChIP motif discovery
    bedtools getfasta \
        -fi "${REF_GENOME_FA}" \
        -bed "${HOMER_OUT}/${sample}/${sample}.bed" \
        -fo "${MEME_OUT}/${sample}/${sample}.fa"

    meme-chip \
        -oc "${MEME_OUT}/${sample}" \
        "${MEME_OUT}/${sample}/${sample}.fa"

    # 8) Intersection with gene annotations
    bedtools intersect \
        -a "${HOMER_OUT}/${sample}/${sample}.bed" \
        -b "${REF_GTF}" \
        > "${HOMER_OUT}/${sample}/${sample}_withAnnotations.bed"
done

module purge
