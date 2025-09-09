# MeRIP-seq Processing Pipeline

## Overview

A Singularity container for processing MeRIP-seq (Methylated RNA Immunoprecipitation sequencing) data. This pipeline provides a reproducible environment for analyzing RNA methylation data, including quality control, read alignment, peak calling, and differential methylation analysis.

## Features

- Quality control with FastQC
- Adapter trimming with Trim Galore
- Read alignment with STAR
- Peak calling with exomePeak2
- Differential methylation analysis
- Support for single-end and paired-end data
- Reproducible environment via Singularity container

## Requirements

- Linux system
- Singularity (version 3.0 or higher)
- Sufficient disk space for data processing
- Minimum 16GB RAM (32GB recommended for large datasets)

## Installation

### Needed software
FastQC (v0.11.9)

Trim Galore (v0.6.7)

STAR (v2.7.10a)

SAMtools (v1.15)

Deeptools (v2.30.0)

MEME Suite (v.5.5.8)

exomePeak2 (v3.16)

R (v4.1.2) with essential packages

## Detailed Usage

### Input Data Preparation

Organize your input files as follows:
```
/data/
  ├── raw_reads/
  │   ├── sample1_IP.fastq.gz
  │   ├── sample1_input.fastq.gz
  │   ├── sample2_IP.fastq.gz
  │   └── sample2_input.fastq.gz
  └── reference/
      ├── hg38/
      │   ├── genome.fa
      │   └── annotation.gtf
      └── STAR_index/
```

## Pipeline Steps

###  1. Raw Data Quality Check(QC)
You can perform quality check on the raw data to assess the sequencing quality.
```
mkdir FastQC
fastqc -o FastQC --noextract -f fastq ./rawdata/*.fq.gz -t 16 >FastQC/fastqc.log 2>&1 &
```
### 2. Adapter Trimming
Based on the quality control results, you can perform appropriate trimming on the raw data.
```
mkdir trim_data
for file in rawdata/*_1.fq.gz; do
    filename=$(echo "$file" | sed -E 's|rawdata/(.*)_1.fq.gz|\1|')
    echo ${filename}
    trim_galore -q 25 --cores 16 --phred33 --fastqc --length 36 -e 0.1 --stringency 3 \
    --paired ./rawdata/${filename}_1.fq.gz ./rawdata/${filename}_2.fq.gz \
    -o ./trim_data >trim_data/trim-${filename}.log 2>&1
done
```
### 3. Read Alignment
Align the quality-controlled FASTQ files to the reference genome, filter for high-quality aligned sequences, and perform sorting.

```
mkdir StarResult
for file in rawdata/*_1.fq.gz; do
    filename=$(echo "$file" | sed -E 's|rawdata/(.*)_1.fq.gz|\1|')
    echo ${filename}
    STAR --runMode alignReads --genomeDir /data/reference/STAR_index \
    --outFileNamePrefix ./StarResult/${filename} \
    --outSAMattributes All  --outSAMtype BAM SortedByCoordinate \
    --readFilesIn ./trim_data/${filename}_1_val_1.fq.gz ./trim_data/${filename}_2_val_2.fq.gz \
    --runThreadN 20 --readFilesCommand zcat \
    --outFilterMultimapNmax 1 --outFilterScoreMinOverLread 0.1 --outFilterMatchNminOverLread 0.1  >4.StarResult/Star_${filename}.log 2>&1
done
```
### 4.  Remove Duplicates
Remove duplicate sequences.
```
samtools fixmate -r -m ./StarResult/${filename}.sort.bam -| samtools sort -@ 20 - | samtools markdup -r -s - ./StarResult/${filename}.no_dup.bam
samtools flagstat ./StarResult/${filename}.no_dup.bam > ./StarResult/${filename}.flagstat.txt
samtools index ./StarResult/${filename}.no_dup.bam
```
### 5. Visualization on Genome Browser
You can visualize these results in the UCSC genome browser or IGV genome browser.
```
mkdir BWFiles
bamCoverage -b ./StarResult/${filename}.no_dup.bam -o BWFiles/${fileName}_forward.bw --filterRNAstrand reverse -p 20 --binSize 1 --normalizeUsing CPM 
bamCoverage -b ./StarResult/${filename}.no_dup.bam -o BWFiles/${fileName}_revers.bw --filterRNAstrand forward -p 20 --binSize 1 --normalizeUsing CPM 
```
### 6. Peak Calling
```
# Download exomePeak2
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rsamtools", "GenomicAlignments", "GenomicRanges", 
                       "GenomicFeatures", "DESeq2", "ggplot2", "mclust", "BSgenome", 
                       "Biostrings", "GenomeInfoDb", "BiocParallel", "IRanges", 
                       "S4Vectors", "rtracklayer", "methods", "stats", 
                       "utils", "BiocGenerics", "magrittr", "speedglm", "splines"))

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("ZW-xjtlu/exomePeak2", build_vignettes = TRUE)

library(exomePeak2)
#Specify File Directories
GENE_ANNO_GTF = "gencode.v38.annotation.gtf"
IP_BAM = c("sample1_IP.bam","sample2_IP.bam")
INPUT_BAM = c("sample1_Input.bam","sample2_Input.bam")
# Peak Calling
sep <- exomePeak2(bam_ip = IP_BAM,
bam_input = INPUT_BAM,
gff_dir = GENE_ANNO_GTF,
genome = "hg38",
paired_end = FALSE)
sep
# You will get the bed files of m6A peaks region
```
### 7. Motif analysis
```
# Get fasta files of m6A peaks region
bed2fasta -s -both -o m6A_peaks_region.fa m6A_peaks_region.bed hg38.fa
# Find motifs in m6A peaks region
streme --oc ./peaks_motif --rna --minw 5 --maxw 5 --thresh 0.05 --align center --p m6A_peaks_region.fa
