

set -euo pipefail


SAMPLE_LIST="./in/ids_RNASeq.txt"

PROJECT_DIR="ppath_rnaseq"
FASTQ_DIR="${PROJECT_DIR}/in/FASTQ"
TRIM_DIR="${PROJECT_DIR}/out/trimmed"
SUMMARY_DIR="${PROJECT_DIR}/out/summary"
BAM_DIR="${PROJECT_DIR}/out/bams"

INDEX_PREFIX="${PROJECT_DIR}/in/hg38/index_files/hg38"

cd "${PROJECT_DIR}"

mkdir -p "${TRIM_DIR}" "${SUMMARY_DIR}" "${BAM_DIR}"

while read sample; do
    mkdir -p "${BAM_DIR}/${sample}"
    # 1) Trim adapters
    ml purge
    ml load GCCcore/8.3.0 Trim_Galore/0.6.5-Java-11.0.2-Python-3.7.4
    trim_galore \
        -j 4 \
        -q 20 \
        --gzip \
        --length 30 \
        --trim-n \
        -o "${TRIM_DIR}" \
        --paired "${FASTQ_DIR}/${sample}_1.fastq.gz" "${FASTQ_DIR}/${sample}_2.fastq.gz"

    # 2) Align & post-process
    ml purge
    ml load HISAT2/2.2.1 SAMtools/1.15
    hisat2 \
        -x "${INDEX_PREFIX}" \
        -1 "${TRIM_DIR}/${sample}_1_val_1.fq.gz" \
        -2 "${TRIM_DIR}/${sample}_2_val_2.fq.gz" \
        -k 1 \
        -p ${SLURM_CPUS_PER_TASK} \
        --no-mixed \
        --no-discordant \
        --summary-file "${SUMMARY_DIR}/${sample}_Summary.txt" \
    | samtools view -@ ${SLURM_CPUS_PER_TASK} -b -F 4 - \
    | samtools fixmate -@ ${SLURM_CPUS_PER_TASK} -O SAM -m - - \
    | samtools sort -@ ${SLURM_CPUS_PER_TASK} -m 5G -O BAM - \
    | samtools markdup -@ ${SLURM_CPUS_PER_TASK} -r - - -O SAM \
    | samtools view -b - \
    > "${BAM_DIR}/${sample}/${sample}.bam"

    ./featureCounts/bin/featureCounts \
        -T ${SLURM_CPUS_PER_TASK} \
        -O \
        -t exon \
        -g gene_id \
        -p \
        -M \
        -a ./in/hg38/hg38.refGene.gtf \
        -o ./out/counts/${sample}/${sample}_counts.txt \
        ./out/bams/${sample}/${sample}.bam

    samtools index -@ ${SLURM_CPUS_PER_TASK} "${BAM_DIR}/${sample}/${sample}.bam"
done < "$SAMPLE_LIST"

