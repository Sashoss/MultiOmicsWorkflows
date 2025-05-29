# ATAC-Seq Processing Pipeline

## Directory Structure

```plaintext

project/
│
├── in/
│   └── references/
│       └── hg38/
│           ├── genome.fa              # hg38 FASTA
│           ├── chrom.sizes            # Chromosome sizes
│           ├── hg38.refGene.gtf       # GTF annotation
│           └── hg38-blacklist.v2.bed  # ENCODE blacklist
│
├── out/
│   ├── bams/                          # Raw alignment BAMs (input)
│   ├── mito_output/                   # chrM metrics & filtered BAMs
│   ├── duplicate_output/              # Duplicate-marked / filtered BAMs
│   ├── blacklist_output/              # Blacklist-filtered BAMs
│   ├── shifted_output/                # ATAC-shifted BAMs
│   ├── visual_output/                 # Normalized bigWig tracks (10 bp bins)
│   ├── macs2_output/                  # MACS2 peak calls
│   ├── foldchange_output/             # Fold-enrichment & p-value tracks
│   ├── frip_output/                   # FRiP statistics
│   ├── homer_output/                  # Peak annotation & HOMER motifs
│   └── meme_output/                   # MEME-ChIP de-novo motifs
│
└── src/
    └── featureCounts/bin/featureCounts # Subread featureCounts binary
```

## <b><u>STEPS</u></b>
WE will skip alignment step here, as its shown in detail in [RNASeq readme file](Notebook/RNASeq/Readme.md). Below are the steps after read alignment. 

### Step 0. chrM Metrics & Removal using [SAMtools](https://www.htslib.org/)
Quantify mitochondrial reads (quality-control proxy) and remove them to avoid inflation of signal-to-noise in downstream analyses.

```bash
samtools idxstats -@ 12 SAMPLE.bam \
  > ./out/mito_output/SAMPLE/SAMPLE.idxstats          # chrM read count

samtools flagstat -@ 12 SAMPLE.bam \
  > ./out/mito_output/SAMPLE/SAMPLE.flagstat          # overall alignment stats

samtools view -@ 12 -h SAMPLE.bam \
  | grep -v '^.*chrM.*$' \
  | samtools sort -@ 12 -O BAM \
    -o ./out/mito_output/SAMPLE/SAMPLE_noChrM.bam

samtools index -@ 12 ./out/mito_output/SAMPLE/SAMPLE_noChrM.bam

```


```plaintext
| Parameter      | Description                            |
| -------------- | -------------------------------------- |
| `idxstats`     | Per-contig read counts (quick QC).     |
| `flagstat`     | Global alignment summary.              |
| `grep -v chrM` | Removes mitochondrial reads in-stream. |
| `-@ 12`        | Use 12 CPU threads.                    |
```


### Step 1. Duplicate Marking & Filtering with [Picard MarkDuplicates](https://broadinstitute.github.io/picard) + SAMtools
PCR duplicates are marked (not removed) first, then high-quality, properly paired reads are retained.

```bash
picard MarkDuplicates \
  INPUT=./out/mito_output/SAMPLE/SAMPLE_noChrM.bam \
  OUTPUT=./out/duplicate_output/SAMPLE/SAMPLE_marked.bam \
  METRICS_FILE=./out/duplicate_output/SAMPLE/SAMPLE.dup.metrics \
  REMOVE_DUPLICATES=false \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=LENIENT

samtools view -@ 12 -h -b -f 2 -F 1548 -q 30 \
  ./out/duplicate_output/SAMPLE/SAMPLE_marked.bam \
| samtools sort -@ 12 -O BAM \
  -o ./out/duplicate_output/SAMPLE/SAMPLE_filtered.bam

samtools index -@ 12 ./out/duplicate_output/SAMPLE/SAMPLE_filtered.bam
```

```plaintext
| Parameter | Description                                                                          |
| --------- | ------------------------------------------------------------------------------------ |
| `-f 2`    | Keep **properly paired** reads only.                                                 |
| `-F 1548` | Remove unmapped/secondary/supplementary/duplicate reads using flag filter mask 1548. |
| `-q 30`   | MAPQ ≥ 30 (high-confidence alignments).                                              |
```

### Step 2. Blacklist Filtering with [BEDTools](https://bedtools.readthedocs.io/en/latest/) intersect
Reads overlapping [ENCODE blacklist regions](https://github.com/igordot/reference-genomes/blob/master/hg38/blacklist.v2.bed) are removed to suppress artefactual signal.

```bash
bedtools intersect -v \
  -abam ./out/duplicate_output/SAMPLE/SAMPLE_filtered.bam \
  -b ./in/references/hg38/hg38-blacklist.v2.bed \
  > ./out/blacklist_output/SAMPLE/SAMPLE_filtered_blacklist.bam

samtools sort -@ 12 -O BAM \
  -o ./out/blacklist_output/SAMPLE/SAMPLE_filtered_blacklist_sorted.bam \
  ./out/blacklist_output/SAMPLE/SAMPLE_filtered_blacklist_sorted.bam

samtools index -@ 12 ./out/blacklist_output/SAMPLE/SAMPLE_filtered_blacklist_sorted.bam
```

```plaintext
| Flag | Description                                                               |
| ---- | ------------------------------------------------------------------------- |
| `-v` | Report ***only*** reads **NOT** overlapping `-b` BED (inverse intersect). |
```

### Step 3. ATAC-shift & Normalized Coverage Tracks with deepTools, alignmentSieve, & bamCoverage
In ATAC-seq the Tn5 transposase binds as a dimer and inserts adapters after cutting a 9-bp staggered duplication. Consequently each read start is displaced **+4 bp on the + strand** and **-5 bp on the – strand** from the true protein–DNA contact site. If this offset is not corrected, aggregate footprints, motif searches, and nucleosome phasing analyses will appear “blurred” and shifted. `alignmentSieve --ATACshift` moves every read by the appropriate amount so that the centre of the transposition event (and thus open chromatin) is accurately represented. Corrects Tn5 offset (±4 bp) and generates 10 bp bin bigWig tracks normalized to BPM (reads per million).

```bash
alignmentSieve \
  --numberOfProcessors $SLURM_CPUS_PER_TASK \
  --ATACshift \
  --bam ./out/blacklist_output/SAMPLE/SAMPLE_filtered_blacklist_sorted.bam \
  -o ./out/shifted_output/SAMPLE/SAMPLE_shifted.bam
| 

samtools sort -@ 12 -O BAM \
            -o ./out/shifted_output/SAMPLE/SAMPLE_shifted_sorted.bam \
            ./out/shifted_output/SAMPLE/SAMPLE_shifted.bam

samtools index -@ 12 ./out/shifted_output/SAMPLE/SAMPLE_shifted_sorted.bam
```

After shifting, we need a quantitative signal that can be visualised in IGV/UCSC or used by downstream QC tools. bamCoverage converts the shifted BAM into a continuous bigWig by counting fragments in fixed windows (here 10 bp). Here, raw coverage is confounded by sequencing depth. The option --normalizeUsing BPM (Bins Per Million mapped reads) divides every bin count by total mapped fragments, making tracks from different samples directly comparable without altering the absolute dynamic range.
In addition, --effectiveGenomeSize excludes Ns and unplaced contigs so that BPM represents reads per mappable million bases, mirroring MACS2’s genome size setting.

```bash
bamCoverage \
  --numberOfProcessors $SLURM_CPUS_PER_TASK \
  --bam ./out/shifted_output/SAMPLE/SAMPLE_shifted.bam \
  --binSize 10 \
  --normalizeUsing BPM \
  --effectiveGenomeSize 2862010578 \
  -o ./out/visual_output/SAMPLE/SAMPLE_shifted.bw

```

This step turns a list of read alignments into an accurate, depth-normalised, visually friendly representation of chromatin accessibility that underpins all subsequent quality checks, inter-sample comparisons, and figure generation.


### Step 4. Peak Calling with [MACS2](https://github.com/macs3-project/MACS) (Broad Mode)
Identify accessible chromatin regions.

```bash
macs2 callpeak \
  -t ./out/shifted_output/SAMPLE/SAMPLE_shifted_sorted.bam \
  -f BAMPE \
  --nomodel \
  --shift -37 \
  --extsize 73 \
  -g 2.86e9 \
  -B --broad --keep-dup all --cutoff-analysis \
  -n SAMPLE \
  --outdir ./out/macs2_output/SAMPLE
```

```plaintext
| Parameter           | Description                                             |
| ------------------- | ------------------------------------------------------- |
| `-f BAMPE`          | Paired-end BAM.                                         |
| `--shift/-extsize`  | Empirical fragment model (-37 / 73 bp) for cut sites.   |
| `--broad`           | Calls **broad** peaks (enhancers, heterochromatin).     |
| `--cutoff-analysis` | Generates peak score distribution for cutoff selection. |
```

### Step 5. Fold-Enrichment & p-Value Tracks
Compute signal-to-noise tracks (FE) and Poisson p-value tracks, then convert to bigWig.

```bash
# Fold-enrichment
macs2 bdgcmp -m FE \
  -t ./out/macs2_output/SAMPLE/SAMPLE_treat_pileup.bdg \
  -c ./out/macs2_output/SAMPLE/SAMPLE_control_lambda.bdg \
  -o ./out/foldchange_output/SAMPLE/SAMPLE_FE.bdg

sort -k1,1 -k2,2n ./out/foldchange_output/SAMPLE/SAMPLE_FE.bdg \
  > ./out/foldchange_output/SAMPLE/SAMPLE_FE.sorted.bdg

bedGraphToBigWig \
  ./out/foldchange_output/SAMPLE/SAMPLE_FE.sorted.bdg \
  ./in/references/hg38/chrom.sizes \
  ./out/foldchange_output/SAMPLE/SAMPLE_FE.bw
```

```plaintext
| File suffix | Meaning                                    |
| ----------- | ------------------------------------------ |
| `_FE.bw`    | Fold-enrichment bigWig.                    |
| `_ppois.bw` | -log<sub>10</sub>(Poisson p-value) bigWig. |
```

### Step 6. FRiP (Fraction of Reads in Peaks) with [featureCounts](https://subread.sourceforge.net/)
Counts reads overlapping called peaks and computes FRiP, a key QC metric.

```bash
awk 'BEGIN{FS=OFS="\t"} NR==1{print "GeneID","Chr","Start","End","Strand"; next} \
     {print $4,$1,$2+1,$3,"."}' \
     ./out/macs2_output/SAMPLE/SAMPLE_broad_peaks.broadPeak \
  > ./out/frip_output/SAMPLE/SAMPLE.saf

featureCounts -p -F SAF \
  -a ./out/frip_output/SAMPLE/SAMPLE.saf \
  -o ./out/frip_output/SAMPLE/SAMPLE-readCountInPeaks.txt \
  ./out/shifted_output/SAMPLE/SAMPLE_shifted.bam
```

```plaintext
| Flag/Option | Description                          |
| ----------- | ------------------------------------ |
| `-p`        | Count **properly paired** fragments. |
| `-F SAF`    | Annotation file is in SAF format.    |
```


### Step 7. Peak Annotation & Known Motifs with [HOMER](http://homer.ucsd.edu/homer/motif/)
Annotate peaks relative to genomic features and searches for enriched known motifs.

```bash
cut -f1-3 ./out/macs2_output/SAMPLE/SAMPLE_broad_peaks.broadPeak \
  > ./out/homer_output/SAMPLE/SAMPLE.bed

annotatePeaks.pl \
  ./out/homer_output/SAMPLE/SAMPLE.bed hg38 \
  > ./out/homer_output/SAMPLE/SAMPLE_annotatedPeaks.txt

findMotifsGenome.pl \
  ./out/homer_output/SAMPLE/SAMPLE.bed hg38 \
  ./out/homer_output/SAMPLE/motifs \
  -preparsedDir ./out/homer_output/SAMPLE/preparsed
```

```plaintext
| Output                 | Purpose                                 |
| ---------------------- | --------------------------------------- |
| `*_annotatedPeaks.txt` | TSS distance, gene mapping, GC content. |
| `motifs/*`             | Known motif enrichment reports.         |
```

### Step 8. De-Novo Motif Discovery with [MEME-ChIP](https://meme-suite.org/meme/tools/meme-chip)
Extracts peak sequences and performs unbiased motif discovery.

```bash
bedtools getfasta \
  -fi ./in/references/hg38/genome.fa \
  -bed ./out/homer_output/SAMPLE/SAMPLE.bed \
  -fo ./out/meme_output/SAMPLE/SAMPLE.fa

meme-chip -oc ./out/meme_output/SAMPLE \
  ./out/meme_output/SAMPLE/SAMPLE.fa
```

```plaintext
| Output directory      | Contents                                   |
| --------------------- | ------------------------------------------ |
| `meme_output/SAMPLE/` | MEME, DREME, CentriMo HTML & text reports. |
```


### Step 9. Peak - Gene Intersection with BEDTools intersect
Overlaps peak BED with gene annotations (GTF converted on-the-fly) to generate custom gene-centric peak files.

```bash
bedtools intersect -a ./out/homer_output/SAMPLE/SAMPLE.bed \
  -b ./in/references/hg38/hg38.refGene.gtf \
  > ./out/homer_output/SAMPLE/SAMPLE_withAnnotations.bed
```

```plaintext
| Result file             | Description                             |
| ----------------------- | --------------------------------------- |
| `*_withAnnotations.bed` | Peaks annotated with overlapping genes. |
```