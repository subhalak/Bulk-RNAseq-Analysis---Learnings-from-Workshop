# Bulk RNAseq Analysis - Learnings from Workshop


## A Complete Reproducible Tutorial / Hands-On Workshop

RNA-seq (RNA sequencing) is a powerful and widely used method to study the transcriptome—the full set of RNA molecules expressed in a given cell type under specific biological conditions. Unlike traditional microarrays, RNA-seq provides a high-resolution, quantitative, and unbiased view of gene expression, enabling the discovery of differentially expressed genes, detection of novel transcripts, alternative splicing, and more.

This repository provides a fully reproducible workflow for Bulk RNA-Seq processing, beginning with raw FASTQ files and ending with a complete count matrix ready for downstream analyses.

## Bulk RNA-seq Processing Pipeline — Guo et al., 2019 (GEO dataset ID: GSE106305)

This tutorial is inspired by the study by Guo et al. (2019), “ONECUT2 is a driver of neuroendocrine prostate cancer,” published in Nature Communications
🔗 https://www.nature.com/articles/s41467-019-11579-6

The study identifies ONECUT2 as a key transcription factor driving neuroendocrine prostate cancer (NEPC)—a highly aggressive and treatment-resistant form of prostate cancer. Integrating transcriptomic and epigenomic datasets, the authors show that ONECUT2 activates neuroendocrine lineage programs, suppresses androgen receptor signaling, and promotes tumor progression. The work establishes ONECUT2 as a potential therapeutic target in NEPC. 

Note: The original analysis integrates transcriptomic profiling across multiple prostate cancer states.

A fully reproducible hands-on workflow for processing raw FASTQ → aligned BAM → gene counts → final count matrix using 20 RNA-seq reads from 8 individual samples (SRR7179504–SRR7179541) from:

<img width="800" height="500" alt="Screenshot 2025-12-10 at 9 34 49 PM" src="https://github.com/user-attachments/assets/4af987ca-654f-468b-b7a1-76c3d1fe4ab1" />

## Dataset Description

Guo et al. utilize multiple publicly available and experimentally generated datasets, including RNA-seq, ChIP-seq, ATAC-seq, and additional prostate cancer transcriptomic datasets from GEO and dbGaP. These datasets include bulk RNA-seq profiles of tumor samples, cell lines, and PDX models.

For this workshop, we demonstrate the workflow using a publicly accessible GEO dataset (Accession ID: [GSE106305](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106305)) containing multiple RNA-seq samples across biological conditions. Raw sequencing data are downloaded via SRA Toolkit, and all files are stored in structured directories for reproducibility.

# Introduction

## Workflow diagram on Bulk RNA seq analysis

Here, we re-process 20 RNA-seq samples publicly available under GSE106305.

The goal: convert raw sequencing data (SRA) → FASTQ → QC → trimmed reads → QC again → Mapping → aligned BAM → gene-level counts → merged count matrix → DESeq2 → GSEA 

<img width="250" height="400" alt="Workflow_diagram" src="https://github.com/user-attachments/assets/185f46cc-3350-47eb-ac63-54bf2e16a938" />

The workflow covers the two major phases of RNA-seq analysis:

## 1. Data preprocessing — converting raw sequencing reads into a count matrix
## 2. Downstream statistical analysis — differential expression, visualization, and biological interpretation

This repository also contains a step-by-step guide for installing required tools, running commands, organizing output files, and performing DESeq2 analysis in R.

## Important Note: Computational requirements to download ~20 FASTQ files?

* Downloading FASTQ files depends on file size, not CPU.

* Typical size per RNA-seq sample:

--Single-end: 1–3 GB

--Paired-end: 3–8 GB

* Public studies vary widely.

--If you download 20 RNA-seq samples, expect:

--Approx required storage - 40–120 GB (realistic range)
--Plan for 150 GB free to be safe.

## 1. Data preprocessing - Methods Overview:

The processing pipeline includes:

* Quality control (FastQC, MultiQC)

* Adapter/quality trimming (Trimmomatic)

* Alignment to reference genome (HISAT2)

* SAM/BAM conversion & indexing (SAMtools)

* Gene-level quantification (featureCounts)

## 2. Downstream statistical analysis:

* Differential Expression Analysis (DESeq2)

* PCA plots, heatmaps, volcano plots

* Pathway & Functional enrichment (ClusterProfiler, GSEA, Reactome, KEGG)

# PART A — Conda Environment Setup

Before running the pipeline, create a dedicated conda environment to ensure consistent versions and reproducible results.

— Install Miniconda (safe, no admin required)

"curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
"bash Miniconda3-latest-MacOSX-arm64.sh"

- Check installation:

"conda --version"

- Create a new environment for the workshop

"conda create -n rnaseq python=3.10 -y"
"conda activate rnaseq"

- Install SRA Toolkit

"conda install -c bioconda sra-tools -y"

## Install all required RNA-seq tools:

"conda install -y -c conda-forge -c bioconda fastqc multiqc trimmomatic hisat2 samtools subread"

## Tools included:

- SRA Toolkit: command-line tools from NCBI for accessing the Sequence Read Archive (SRA) using prefetch, fastq-dump, fasterq-dump etc

- FastQC: Assess quality of raw FASTQ reads

- MultiQC: aggregates and summarizes output from many different analysis tools like FastQC

- Trimmomatic: Trim adapters & low-quality bases

- HISAT2: Splice-aware alignment

- SAMtools: Handle SAM/BAM files

- Subread (featureCounts): Gene-level quantification

# PART B — Creating Directory Structure

## Organize your working directory:

Run the following command as follows

"mkdir rnaseq"
"cd rnaseq"
"mkdir fastq"

# PART C — Downloading GEO Datasets (SRA)

## Install SRA Toolkit:

"sudo apt install sra-toolkit"
"prefetch SRR7179504"

## Download FASTQ files using fasterq-dump:

"fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip SRR7179504.sra"

These are NCBI SRA Toolkit commands for downloading sequencing data.

"cd fastq"

- Download all 20 .sra files

- Convert them into 20 FASTQ files (.fastq.gz) inside your fastq/ folder using "python3 file_name.py"

Note: 

- prefetch SRRxxxx → download files

- fastq-dump → convert .sra to .fastq.gz

## Typical size:

- Each .sra file: ~1–3 GB

- After converting to fastq: ~1–5 GB per sample

## To verify: 

"ls fastq"

we should see files like:

"SRR7179504_1.fastq.gz"
"SRR7179504_2.fastq.gz"

If we see all 20 samples → we are READY for the workshop.

# PART D — Downloading Reference Genome & Annotation Files

## Download HISAT2 GRCh38 index:

"wget -P $REFERENCE https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz"
"tar -xzf grch38_genome.tar.gz"

## Download GTF annotation:

"wget -P $REFERENCE https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz"
"gunzip Homo_sapiens.GRCh38.109.gtf.gz"

##STEP 1 — Quality Control (FastQC + MultiQC)

## Run FastQC:

"fastqc $FASTQ/*_1.fastq $FASTQ/*_2.fastq -o $FASTQC_Results -t 10"

## Generate summary report:

"multiqc $FASTQC_Results -o $FASTQC_Results"

# STEP 2 — Trimming (Trimmomatic)

Example command:

"trimmomatic PE input_1.fastq input_2.fastq out_1_paired.fq out_1_unpaired.fq out_2_paired.fq out_2_unpaired.fq SLIDINGWINDOW:4:20 MINLEN:36"

# STEP 3 — Alignment with HISAT2

## Single sample:

"hisat2 -p 20 -x $REFERENCE/grch38/genome -1 sample_1.fastq -2 sample_2.fastq -S ALIGN/sample_aligned.sam"

## Loop over all samples:

"for fq1 in $FASTQ/*_1.fastq; do
    fq2=${fq1/_1.fastq/_2.fastq}
    base=$(basename $fq1 _1.fastq)
    hisat2 -p 20 -x $REFERENCE/grch38/genome -1 $fq1 -2 $fq2 -S $ALIGN/${base}_aligned.sam
done"

# STEP 4 — SAM → BAM Conversion (SAMtools)

## Convert, sort, and index:

"for sam in $ALIGN/*_aligned.sam; do
    base=$(basename $sam _aligned.sam)
    samtools view -@ 20 -bS $sam | samtools sort -@ 20 -o $BAM/${base}_sorted.bam
    samtools index $BAM/${base}_sorted.bam
done"

# STEP 5 — Read Counting (featureCounts)

"featureCounts -S 2 -a Homo_sapiens.GRCh38.114.gtf -o quants/featurecounts.txt sample.bam"

# Loop for all BAM files:

"for f in $BAM/*_sorted.bam; do
    featureCounts -T 20 -p -a $REFERENCE/Homo_sapiens.GRCh38.109.gtf -o $COUNTS/$(basename ${f%_sorted.bam})_counts.txt "$f"
done"

## FINAL — Merge All Count Files

"mkdir MERGED_COUNTS"

"awk 'NR>1 {print $1}' $COUNTS/SRR11262284_counts.txt > MERGED_COUNTS/geneids.txt"

"for f in $COUNTS/*_counts.txt; do
    sample=$(basename "$f" _counts.txt)
    awk 'NR>1 {print $7}' "$f" > MERGED_COUNTS/${sample}_col.txt
done"

"paste MERGED_COUNTS/geneids.txt MERGED_COUNTS/*_col.txt > MERGED_COUNTS/merged_counts.txt"

# Downstream Analysis — DESeq2 in R

This section covers:

* Importing count matrix

* Running DESeq2

* PCA / heatmap / volcano plot

* GO/KEGG enrichment

# Differential Expression Analysis (DESeq2)

Install required packages:

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

Run - 
# End of Tutorial

This README provides a full, reproducible workflow for bulk RNA-seq—from raw FASTQ files to biological interpretation of differentially expressed genes.
