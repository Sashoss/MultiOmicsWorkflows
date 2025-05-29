# RNA-Seq Processing Pipeline

## Directory Structure

```plaintext
project/
│
├── in/
│   ├── FASTQ/                    # Raw FASTQ files organized by sample folder
│   └── ids_RNASeq.txt            # List of sample IDs, one per line
│
├── out/
│   ├── bams/                     # Aligned BAM files
│   ├── counts/                   # Gene-level counts
│   ├── summary/                  # Summary stats
│   └── trimmed/                  # Trimmed FASTQ files
│
├── scripts/
│   └── count_inline.py           # Script to count lines in piped streams
│
├── Reference/
│   └── hg38/index_files/         # HISAT2 genome index files
│
├── Data/References/
│   └── hg38/hg38.refGene.gtf     # GTF annotation file
│
└── featureCounts/                # featureCounts binary
```

<hr>

## <b><u>STEPS</u></b>

<hr>

### Step 0. Quality Check using <b>[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)</b> 
FastQC provides a comprehensive quality assessment of raw and trimmed reads. This step helps identify adapter contamination, low-quality regions, duplication levels, GC content bias, and more.

```bash
fastqc -t 12 \
  -o ./out/qc/raw/Sample001 \
  ./in/FASTQ/Sample001/Sample001_RNASeq_R1.fastq.gz \
  ./in/FASTQ/Sample001/Sample001_RNASeq_R2.fastq.gz
```

```plaintext
| Parameter | Description                          |
| --------- | ------------------------------------ |
| `-t 12`   | Number of threads to use.            |
| `-o`      | Output directory for FastQC reports. |
```

<hr>

### Step 1. Read trimming using <b>[Trim Galore](https://github.com/FelixKrueger/TrimGalore)</b>
Remove adapter sequences, trim low-quality bases, and filter out short or ambiguous reads to improve alignment efficiency and reduce downstream noise.
- Task
   - Trimming adapters: Removes sequencing adapters that may interfere with mapping.
   - Quality filtering: Bases below Q20 are trimmed from the ends.
   - Length filtering: Discards reads shorter than 30 bp after trimming.
   - Removes Ns: Gets rid of ambiguous base calls (N).

```bash
trim_galore -j 12 -q 20 --gzip --length 30 --trim-n -o ./out/trimmed/Sample001 --paired \
  ./in/FASTQ/Sample001/Sample001_RNASeq_R1.fastq.gz \
  ./in/FASTQ/Sample001/Sample001_RNASeq_R2.fastq.gz
```

```plaintext
| Parameter     | Description                                       |
| ------------- | ------------------------------------------------- |
| `-j 12`       | Number of cores used for parallel processing.     |
| `-q 20`       | Minimum quality score for trimming.               |
| `--length 30` | Minimum length of reads to retain after trimming. |
| `--trim-n`    | Remove Ns at the end of reads.                    |
| `--paired`    | Indicates paired-end read trimming.               |
| `--gzip`      | Compresses the output FASTQ files.                |
| `-o`          | Output directory for trimmed reads.               |
```

<hr>


### Step 2. Download [hg38 human genome](https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/) and create [Hisat2](https://daehwankimlab.github.io/hisat2/) index files

```bash
# Download hg38 genome (FASTA format)
wget https://hgdownload.soe.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

# Create a directory to store index files
mkdir -p ../Reference/hg38/index_files

# Build HISAT2 index
hisat2-build -p 12 hg38.fa ../Reference/hg38/index_files/hg38

```

### Step 3. Read alignment to hg38 human genome using Hisat2
Hisat2 argumentsshown below ensures high-confidence mapping of reads to the genome, eliminating poor or ambiguous mappings. Paired-end reads are treated strictly to reduce noise in downstream analysis.

```bash
hisat2 \
  -x ../Reference/hg38/index_files/hg38 \
  -1 ./out/trimmed/Sample001/Sample001_RNASeq_R1_val_1.fq.gz \
  -2 ./out/trimmed/Sample001/Sample001_RNASeq_R2_val_2.fq.gz \
  -k 1 \
  -p 12 \
  --no-mixed \
  --no-discordant \
  --summary-file ./out/summary/Sample001/Sample001_RNASeq_Summary.txt \
  -S ./out/summary/Sample001/Sample001_RNASeq.sam

```

| Parameter         | Description                                                                  |
| ----------------- | ---------------------------------------------------------------------------- |
| `-x`              | Path to the **HISAT2 index** for the reference genome (e.g., hg38).          |
| `-1`, `-2`        | Paired-end **trimmed FASTQ** files to be aligned.                            |
| `-k 1`            | Reports **only the best alignment** per read (reduces ambiguity).            |
| `-p`              | Number of **threads** to use (faster performance).                           |
| `--no-mixed`      | Disallows **unaligned pairs** where only one read maps (stricter quality).   |
| `--no-discordant` | Disallows **abnormally mapped pairs** (e.g., wrong orientation or distance). |
| `--summary-file`  | Saves a **summary** of alignment statistics for this sample.                 |


<hr>

### Step 4: Convert SAM to BAM and Remove Unmapped Reads using [SAMtools](https://www.htslib.org/)
Filters out unmapped reads (-F 4) and converts to BAM.

```bash
samtools view -@ 12 -b -F 4 \
  ./out/summary/Sample001/Sample001_RNASeq.sam \
  -o ./out/summary/Sample001/Sample001_RNASeq.filtered.bam
```

```plaintext
| Parameter | Description                    |
| --------- | ------------------------------ |
| `-b`      | Output in BAM format.          |
| `-F 4`    | Exclude unmapped reads.        |
| `-@ 12`   | Use 12 threads for processing. |
```

<hr>

### Step 5: Fix Mate Pair Information
Ensures proper pairing information in the BAM file, which is essential for downstream sorting and duplicate removal.

```bash
samtools fixmate -@ 12 -m \
  ./out/summary/Sample001/Sample001_RNASeq.filtered.bam \
  ./out/summary/Sample001/Sample001_RNASeq.fixmate.bam
```

```plaintext
| Parameter | Description                          |
| --------- | ------------------------------------ |
| `-m`      | Marks secondary/supplementary reads. |
| `-@ 12`   | Use 12 threads.                      |
```

<hr>

### Step 6: Sort BAM File
Sorts reads by genomic coordinate. Required for duplicate removal.

```bash
samtools sort -@ 12 -m 5G -O BAM \
  -o ./out/summary/Sample001/Sample001_RNASeq.sorted.bam \
  ./out/summary/Sample001/Sample001_RNASeq.fixmate.bam
```

```plaintext
| Parameter | Description                       |
| --------- | --------------------------------- |
| `-m 5G`   | Maximum memory per thread (5 GB). |
| `-O BAM`  | Output format set to BAM.         |
| `-@ 12`   | Use 12 threads.                   |
```

<hr>

### Step 7: Remove PCR Duplicates
Removes duplicate reads to reduce false quantification.

```bash
samtools markdup -@ 12 -r \
  ./out/summary/Sample001/Sample001_RNASeq.sorted.bam \
  ./out/summary/Sample001/Sample001_RNASeq.dedup.sam \
  -O SAM
```

```plaintext
| Parameter | Description                             |
| --------- | --------------------------------------- |
| `-r`      | Remove (not just mark) duplicate reads. |
| `-@ 12`   | Use 12 threads.                         |
| `-O SAM`  | Output format is SAM.                   |
```

<hr>

### Step 8: Convert SAM to Final BAM file

```bash
samtools view -@ 12 -b \
  ./out/summary/Sample001/Sample001_RNASeq.dedup.sam \
  -o ./out/bams/Sample001/Sample001_RNASeq.bam
```

```plaintext
| Parameter | Description               |
| --------- | ------------------------- |
| `-b`      | Output format set to BAM. |
| `-@ 12`   | Use 12 threads.           |
```

<hr>