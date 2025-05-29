# CHIP-Seq Processing Pipeline

## Directory Structure

```plaintext
project/
│
├── in/
│   └── hg38/
│       ├── genome.fa              # hg38 FASTA
│       ├── hg38.refGene.gtf       # GTF annotation
│       ├── hg38-blacklist.v2.bed  # ENCODE blacklist
│       └── chrom.sizes            # Chromosome sizes
│
├── out/
│   ├── bams/                       # Raw alignment BAMs (input)
│   ├── qc/                         # Pre-peak QC (flagstat, idxstats, insert size)
│   ├── filtered_bams/              # chrM-removed, deduped, MAPQ-filtered BAMs
│   ├── macs2_output/               # MACS2 peak calls
│   ├── tracks/                     # Browser tracks (bigWig)
│   ├── frip_output/                # FRiP statistics
│   ├── homer_output/               # Peak annotation & HOMER motifs
│   ├── meme_output/                # de novo motif discovery
│   └── idr_output/                 # IDR output (if used downstream)
│
└── src/
    └── bedGraphToBigWig           # Custom bedGraph→bigWig converter
```

<hr>
<br>

## <b><u>STEPS</u></b>
WE will skip alignment step here, as its shown in detail in [RNASeq readme file](Notebook/RNASeq/Readme.md). Below are the steps after read alignment. 

### Step 0. Pre-peak QC: flagstat, idxstats, insert-size metrics

- Flagstat: overall alignment summary

```bash
samtools flagstat SRRXXXXXXX.bam \
  > ./out/qc/SRRXXXXXXX/SRRXXXXXXX.flagstat
```

- Idxstats: per-chromosome read counts

```bash
samtools idxstats SRRXXXXXXX.bam \
  > ./out/qc/SRRXXXXXXX/SRRXXXXXXX.idxstats
```

- Insert-size metrics (Picard)

```bash
java -jar picard.jar CollectInsertSizeMetrics \
  INPUT=SRRXXXXXXX.bam \
  OUTPUT=./out/qc/SRRXXXXXXX/SRRXXXXXXX_insert_size_metrics.txt \
  HISTOGRAM_FILE=./out/qc/SRRXXXXXXX/SRRXXXXXXX_insert_size_hist.pdf \
  VALIDATION_STRINGENCY=LENIENT
```


```plaintext
| Metric                     | Purpose                            |
| -------------------------- | ---------------------------------- |
| `flagstat`                 | Global alignment summary           |
| `idxstats`                 | Quick per-chromosome QC            |
| `CollectInsertSizeMetrics` | Fragment length distribution check |
```

<hr>
<br>

### Step 1. Remove chrM, mark duplicates, filter MAPQ 
- Remove mitochondrial reads

```bash
samtools view -h SRR.treatment.bam \
  | grep -v '^.*chrM.*$' \
  > SRR_noChrM.unsorted.bam

samtools sort -O BAM \
  -o ./out/filtered_bams/SAMPLE/SAMPLE_noChrM.bam \
  SRR_noChrM.unsorted.bam

samtools index ./out/filtered_bams/SAMPLE/SAMPLE_noChrM.bam
```

- Mark (and remove) duplicates

```bash
java -jar picard.jar MarkDuplicates \
  INPUT=./out/filtered_bams/SAMPLE/SAMPLE_noChrM.bam \
  OUTPUT=./out/filtered_bams/SAMPLE/SAMPLE_dedup.bam \
  METRICS_FILE=./out/qc/SAMPLE/SAMPLE_dup_metrics.txt \
  REMOVE_DUPLICATES=true \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT
```

- Keep only high-quality, properly reads

```bash
samtools view -@ 10 -b \
    -f 2 \       # if paired
    -q 30 \      
    ./out/filtered_bams/SAMPLE/SAMPLE_dedup.bam \
  | samtools sort -O BAM \
      -o ./out/filtered_bams/SAMPLE/SAMPLE_filtered.bam

samtools index ./out/filtered_bams/SAMPLE/SAMPLE_filtered.bam
```

<hr>
<br>

### Step 2. MACS2 Peak Calling

```bash
macs2 callpeak \
  -t ./out/filtered_bams/SAMPLE/SAMPLE_filtered.bam \
  -c ./out/filtered_bams/SAMPLE/default_input_filtered.bam \
  -f BAM \
  --call-summits \
  -g hs \
  -B --SPMR \              
  --cutoff-analysis \
  -n SAMPLE \
  --outdir ./out/macs2_output/SAMPLE
```

```plaintext
| Parameter           | Description                                   |
| ------------------- | --------------------------------------------- |
| `--call-summits`    | Report highest-resolution summit per peak     |
| `-g hs`             | Human genome size                             |
| `-B --SPMR`         | Output bedGraph normalized to RPM             |
| `--cutoff-analysis` | Generate score distribution for cutoff choice |
```

<hr>
<br>

### Step 3. Lets generate Tracks to visualize in browser (bedGraph → bigWig)

- Sort the bedGraph

```bash
sort -k1,1 -k2,2n \
  ./out/macs2_output/SAMPLE/SAMPLE_treat_pileup.bdg \
  > ./out/tracks/SAMPLE/SAMPLE_pileup.sorted.bdg
```

- Convert to bigWig (for easy visualization)

```bash
./src/bedGraphToBigWig \
  ./out/tracks/SAMPLE/SAMPLE_pileup.sorted.bdg \
  ./in/hg38/chrom.sizes \
  ./out/tracks/SAMPLE/SAMPLE_pileup.bw
```

<hr>
<br>

### Step 4. Filter Peaks Against ENCODE Blacklist

```bash
bedtools intersect -v \
  -a ./out/macs2_output/SAMPLE/SAMPLE_peaks.narrowPeak \
  -b ./in/hg38/hg38-blacklist.v2.bed \
  > ./out/macs2_output/SAMPLE/SAMPLE_noblacklist.narrowPeak
```

Here, -v reports only those peaks not overlapping blacklist regions.

<hr>
<br>

### Step 5. FRiP Score (Fraction of Reads in Peaks)

- Create SAF annotation for peaks

```bash
awk 'BEGIN{FS=OFS="\t"} 
    { if (NF==10) print $4, $1, $2+1, $3, "+"; }' \
  ./out/macs2_output/SAMPLE/SAMPLE_noblacklist.narrowPeak \
  > ./out/frip_output/SAMPLE/SAMPLE.saf
```

- Count reads in peaks

```bash
featureCounts \
  -p \                   # properly paired fragments
  -F SAF \
  -a ./out/frip_output/SAMPLE/SAMPLE.saf \
  -o ./out/frip_output/SAMPLE/SAMPLE_readCountInPeaks.txt \
  ./out/filtered_bams/SAMPLE/SAMPLE_filtered.bam
```

<hr>
<br>

### Step 6. Peak Annotation & Known Motif Enrichment (HOMER)

- Extract peak coordinates

```bash
cut -f1-3 \
  ./out/macs2_output/SAMPLE/SAMPLE_noblacklist.narrowPeak \
  > ./out/homer_output/SAMPLE/SAMPLE.bed
```

- Annotate peaks

```bash
annotatePeaks.pl \
  ./out/homer_output/SAMPLE/SAMPLE.bed \
  hg38 \
  > ./out/homer_output/SAMPLE/SAMPLE_annotatedPeaks.txt
```

- Find known motifs

```bash
findMotifsGenome.pl \
  ./out/homer_output/SAMPLE/SAMPLE.bed \
  hg38 \
  ./out/homer_output/SAMPLE \
  -preparsedDir ./out/homer_output/SAMPLE/preparsed
```

<hr>
<br>

### Step 7. De Novo Motif Discovery (MEME-ChIP)

- Extract peak sequences

```bash
bedtools getfasta \
  -fi ./in/hg38/genome.fa \
  -bed ./out/homer_output/SAMPLE/SAMPLE.bed \
  -fo ./out/meme_output/SAMPLE/SAMPLE.fa
```

- Run MEME-ChIP

```bash
meme-chip \
  -oc ./out/meme_output/SAMPLE \
  ./out/meme_output/SAMPLE/SAMPLE.fa
```

<hr>
<br>

### Step 8. Peak–Gene Intersection
Here we assign peaks to their corresponding genes.

```bash
bedtools intersect \
  -a ./out/homer_output/SAMPLE/SAMPLE.bed \
  -b ./in/hg38/hg38.refGene.gtf \
  > ./out/homer_output/SAMPLE/SAMPLE_withAnnotations.bed
```

<hr>
<br>