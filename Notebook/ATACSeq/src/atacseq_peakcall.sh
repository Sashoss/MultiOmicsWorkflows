#!/bin/bash
#SBATCH --account=gdleelab
#SBATCH --partition=himem               
#SBATCH --mem=100G  
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name atacseq_pipeline
#SBATCH --partition=himem
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --mail-user=ambuj.kumar@nationwidechildrens.org
#SBATCH --mail-type=ALL

module load miniforge3/24.3.0
conda activate bed2bw

set -euo pipefail

samples=(
    "SRR21418681"
    "SRR21418691"
    "SRR21418697"
)

PROJECT_DIR="/home/gdwanglab/axk201/personal_projects/MultiOmicsWorkflow/Notebook/ATACSeq"
BAM_DIR="${PROJECT_DIR}/out/bams"
MITO_OUT="${PROJECT_DIR}/out/mito_output"
DUPLICATE_OUT="${PROJECT_DIR}/out/duplicate_output"
BLACKLIST_OUT="${PROJECT_DIR}/out/blacklist_output"
SHIFTED_OUT="${PROJECT_DIR}/out/shifted_output"
VISUAL_OUT="${PROJECT_DIR}/out/visual_output"
MACS2_OUT="${PROJECT_DIR}/out/macs2_output"
FOLDCHANGE_OUT="${PROJECT_DIR}/out/foldchange_output"
FRIP_OUT="${PROJECT_DIR}/out/frip_output"
HOMER_OUT="${PROJECT_DIR}/out/homer_output"
MEME_OUT="${PROJECT_DIR}/out/meme_output"

REF_GENOME_FA="${PROJECT_DIR}/in/hg38/genome.fa"
REF_GTF="${PROJECT_DIR}/in/hg38/hg38.refGene.gtf"
REF_chr_sizes="${PROJECT_DIR}/in/hg38/chrom.sizes"
BLACKLIST_BED="${PROJECT_DIR}/in/hg38/hg38-blacklist.v2.bed"

cd "${PROJECT_DIR}"

for sample in "${samples[@]}"; do
    BAM="${BAM_DIR}/${sample}/${sample}.bam"

    # check treatment BAM
    if [[ ! -f "${BAM}" ]]; then
        echo "[WARNING] Missing BAM for sample ${sample}: ${BAM}"
        continue
    fi

    # make per‐sample dirs
    for d in \
        "${MITO_OUT}" "${DUPLICATE_OUT}" "${BLACKLIST_OUT}" \
        "${SHIFTED_OUT}" "${VISUAL_OUT}" "${MACS2_OUT}" \
        "${FOLDCHANGE_OUT}" "${FRIP_OUT}" "${HOMER_OUT}" "${MEME_OUT}"
    do
        mkdir -p "${d}/${sample}"
    done

    # 0) ChrM metrics and remove chrM
    module load SAMtools/1.15
    samtools idxstats -@ 12 "${BAM}" > "${MITO_OUT}/${sample}/${sample}.idxstats"
    grep "^chrM" "${MITO_OUT}/${sample}/${sample}.idxstats" \
         > "${MITO_OUT}/${sample}/${sample}_chrM.log"
    samtools flagstat -@ 12 "${BAM}" \
         > "${MITO_OUT}/${sample}/${sample}.flagstat"
    samtools view -@ 12 -h "${BAM}" \
         | grep -v "^.*chrM.*$" \
         | samtools sort -O BAM \
                        -o "${MITO_OUT}/${sample}/${sample}_noChrM.bam"
    samtools index -@ 12 "${MITO_OUT}/${sample}/${sample}_noChrM.bam"

    # 1) Mark duplicates & filter
    ml load picard/2.21.6-Java-11.0.2
    java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
        QUIET=true \
        INPUT="${MITO_OUT}/${sample}/${sample}_noChrM.bam" \
        OUTPUT="${DUPLICATE_OUT}/${sample}/${sample}_marked.bam" \
        METRICS_FILE="${DUPLICATE_OUT}/${sample}/${sample}.dup.metrics" \
        REMOVE_DUPLICATES=false \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=LENIENT \
        TMP_DIR="."

    module load SAMtools/1.15
    samtools view -@ 12 -h -b -F 1548 -q 30 \
        "${DUPLICATE_OUT}/${sample}/${sample}_marked.bam" \
      | samtools sort -@ 12 -O BAM \
                      -o "${DUPLICATE_OUT}/${sample}/${sample}_filtered.bam"
    samtools index -@ 12 "${DUPLICATE_OUT}/${sample}/${sample}_filtered.bam"

    # 2) Remove ENCODE blacklist regions
    module purge
    ml load GCC/10.3.0
    ml load BEDTools/2.30.0
    bedtools intersect -v \
        -abam "${DUPLICATE_OUT}/${sample}/${sample}_filtered.bam" \
        -b "${BLACKLIST_BED}" \
      > "${BLACKLIST_OUT}/${sample}/${sample}_filtered_blacklist.bam"

    module load SAMtools/1.15
    samtools sort -@ 12 -O BAM \
        -o "${BLACKLIST_OUT}/${sample}/${sample}_filtered_blacklist_sorted.bam" \
        "${BLACKLIST_OUT}/${sample}/${sample}_filtered_blacklist.bam"
    samtools index -@ 12 "${BLACKLIST_OUT}/${sample}/${sample}_filtered_blacklist_sorted.bam"

    # 3) ATAC‐shift & bigWig
    module purge
    ml load GCC/9.3.0
    
    ml load OpenMPI/4.0.3
    ml load deepTools/3.3.1-Python-3.8.2
    
    alignmentSieve \
        --numberOfProcessors $SLURM_CPUS_PER_TASK \
        --ATACshift \
        --blackListFileName "${BLACKLIST_BED}" \
        --bam "${BLACKLIST_OUT}/${sample}/${sample}_filtered_blacklist_sorted.bam" \
        -o ${SHIFTED_OUT}/${sample}/${sample}_filtered_blacklist_shifted.bam
    samtools index -@ 12 "${SHIFTED_OUT}/${sample}/${sample}_filtered_blacklist_shifted.bam"

    module load SAMtools/1.15
    samtools sort -@ 12 -O BAM \
            -o "${SHIFTED_OUT}/${sample}/${sample}_filtered_blacklist_shifted_sorted.bam" \
            "${SHIFTED_OUT}/${sample}/${sample}_filtered_blacklist_shifted_sorted.bam"
    samtools index -@ 12 "${SHIFTED_OUT}/${sample}/${sample}_filtered_blacklist_shifted_sorted.bam"


    bamCoverage \
        --numberOfProcessors ${SLURM_CPUS_PER_TASK} \
        --binSize 10 \
        --normalizeUsing BPM \
        --effectiveGenomeSize 2862010578 \
        --bam "${SHIFTED_OUT}/${sample}/${sample}_filtered_blacklist_shifted_sorted.bam" \
        -o "${VISUAL_OUT}/${sample}/${sample}_filtered_blacklist_shifted.bw"

    # 4) MACS2 broad peak calling
    module purge
    ml load GCC/10.3.0
    ml load OpenMPI/4.1.1
    ml load MACS2/2.2.7.1
    macs2 callpeak \
        -t "${SHIFTED_OUT}/${sample}/${sample}_filtered_blacklist_shifted_sorted.bam" \
        -f BAMPE \
        --nomodel \
        --shift -37 \
        --extsize 73 \
        -g 2862010578 \
        -B --broad \
        --keep-dup all \
        --cutoff-analysis \
        -n "${sample}" \
        --outdir "${MACS2_OUT}/${sample}"

    # 5) Fold‐enrichment & p‐value bigWigs
    macs2 bdgcmp \
        -t "${MACS2_OUT}/${sample}/${sample}_treat_pileup.bdg" \
        -c "${MACS2_OUT}/${sample}/${sample}_control_lambda.bdg" \
        -m FE \
        -o "${FOLDCHANGE_OUT}/${sample}/${sample}_FE.bdg"

    sort -k1,1 -k2,2n "${FOLDCHANGE_OUT}/${sample}/${sample}_FE.bdg" \
         > "${FOLDCHANGE_OUT}/${sample}/${sample}_FE.sorted.bdg"
         
    bedGraphToBigWig \
        "${FOLDCHANGE_OUT}/${sample}/${sample}_FE.sorted.bdg" \
        "${REF_chr_sizes}" \
        "${FOLDCHANGE_OUT}/${sample}/${sample}_FE.sorted.bw"

    macs2 bdgcmp \
        -t "${MACS2_OUT}/${sample}/${sample}_treat_pileup.bdg" \
        -c "${MACS2_OUT}/${sample}/${sample}_control_lambda.bdg" \
        -m ppois \
        -o "${FOLDCHANGE_OUT}/${sample}/${sample}_ppois.bdg"

    sort -k1,1 -k2,2n "${FOLDCHANGE_OUT}/${sample}/${sample}_ppois.bdg" \
         > "${FOLDCHANGE_OUT}/${sample}/${sample}_ppois.sorted.bdg"
         
    bedGraphToBigWig \
        "${FOLDCHANGE_OUT}/${sample}/${sample}_ppois.sorted.bdg" \
        "${REF_chr_sizes}" \
        "${FOLDCHANGE_OUT}/${sample}/${sample}_ppois.sorted.bw"

    # 6) FRiP score
    awk 'BEGIN{FS=OFS="\t"} NR==1{print "GeneID","Chr","Start","End","Strand"; next} \
        {print $4, $1, $2+1, $3, "."}' \
        "${MACS2_OUT}/${sample}/${sample}_peaks.broadPeak" \
      > "${FRIP_OUT}/${sample}/${sample}.saf"

    ./src/featureCounts/bin/featureCounts \
        -p \
        -F SAF \
        -a "${FRIP_OUT}/${sample}/${sample}.saf" \
        -o "${FRIP_OUT}/${sample}/${sample}_readCountInPeaks.txt" \
        "${SHIFTED_OUT}/${sample}/${sample}_filtered_blacklist_shifted.bam"

    # 7) HOMER annotation & motif finding
    module purge
    module load homer/4.11.1
    cut -f1-3 \
        "${MACS2_OUT}/${sample}/${sample}_peaks.broadPeak" \
      > "${HOMER_OUT}/${sample}/${sample}.bed"

    annotatePeaks.pl \
        "${HOMER_OUT}/${sample}/${sample}.bed" \
        hg38 \
      > "${HOMER_OUT}/${sample}/${sample}_annotatedPeaks.txt"

    findMotifsGenome.pl \
        "${HOMER_OUT}/${sample}/${sample}.bed" \
        hg38 \
        "${HOMER_OUT}/${sample}/motifs" \
        -preparsedDir "${HOMER_OUT}/${sample}/preparsed"

    # 8) MEME‐ChIP
    module purge
    ml load GCC/10.3.0
    module load BEDTools/2.30.0
    bedtools getfasta \
        -fi "${REF_GENOME_FA}" \
        -bed "${HOMER_OUT}/${sample}/${sample}.bed" \
        -fo "${MEME_OUT}/${sample}/${sample}.fa"

    module purge
    ml load GCC/9.3.0
    ml load MEME-suite/5.5.2-GCC-9.3.0
    meme-chip \
        -meme-p 12 \
        -oc "${MEME_OUT}/${sample}" \
        "${MEME_OUT}/${sample}/${sample}.fa"

    # 9) Intersection with gene annotations
    module purge
    ml load GCC/10.3.0
    module load BEDTools/2.30.0
    bedtools intersect \
        -a "${HOMER_OUT}/${sample}/${sample}.bed" \
        -b "${REF_GTF}" \
        > "${HOMER_OUT}/${sample}/${sample}_withAnnotations.bed"
done



