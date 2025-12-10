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

## Install Miniconda (safe, no admin required)

"curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh"
"bash Miniconda3-latest-MacOSX-arm64.sh"

## Check installation:

"conda --version"

## Create a new environment for the workshop

"conda create -n rnaseq python=3.10 -y"
"conda activate rnaseq"

## Install SRA Toolkit

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

we will see this verbose after running fastq dump command

<img width="609" height="98" alt="Screenshot 2025-12-10 at 11 37 47 PM" src="https://github.com/user-attachments/assets/24849d64-ebc4-4a43-963b-2eea680c455a" />


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

## Quality Control (FastQC + MultiQC)

## Step 1: Run FastQC:

"fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8"

This will give u html files for each of ur fastq.gz but we need one file for all samples, n for that we do multiqc

## Generate summary report:

"multiqc fastqc_results/ -o multiqc_report/" 

For trimming, we use trimmomatic (installation of which is challenging - try installing it by ur own - hint (githubs links) and run the below command just for one sample

## STEP 2 — Trimming (Trimmomatic)

"java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 fastq/SRR7179504.fastq.gz fastq/SRR7179504_trimmed.fastq.gz TRAILING:10 -phred33"   

After trimming, run FastQC again to compare metrics.

# PART D — Downloading Reference Genome & Annotation Files

## STEP 3 — Download HISAT2 & Alignment with HISAT2 using GRCh38 index:

"wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz"
"tar -xvzf grch38_genome.tar.gz"

"sudo apt install hisat2"
"sudo apt install samtools"

"hisat2 -q -x grch38/genome -U fastq/sample1.fastq.gz | samtools sort -o sample1.bam | samtools index sample1.bam"

## try this (will take like 20-25 mins)  ## single sample command

Now Make a .sh script of the code given below and run ./script_name.sh (this will automate all the alignment process)

# STEP 5 — Read Counting (featureCounts)

As a result, we will see 8 .bam and 8 corresponding .bai files

Homo_sapiens.GRCh38.99.gtf.gz ← given in the tutorial; u need to download the latest version

featureCounts -S 2 -a Homo_sapiens.GRCh38.114.gtf -o quants/featurecounts.txt sample.bam  ## to get counts for from .bam file

# Loop for all BAM files:

Run this ./script_name.sh file to run all 8 pairs of files to generate featurecounts

./script_name.sh  (this will automate estimation of counts by featurecounts taking each .bam and .bai file)

## FINAL — Merge All Count Files

#!/usr/bin/env python
#coding: utf-8

#In[1]:

import os
import glob
import pandas as pd
import time

path = "path_to_folder_where_featurecount_output_is_stored"
files = glob.glob(os.path.join(path, "*.txt"))

print("Files found:", files)

all_counts = []

for file in files:
    start_time = time.time()
    df = pd.read_csv(file, sep="\t", comment="#")
    
    sample_name = os.path.basename(file).replace("_featurecounts.txt", "")
    df = df[["Geneid", df.columns[-1]]]
    df.rename(columns={df.columns[-1]: sample_name}, inplace=True)
    
    all_counts.append(df)
    
    elapsed = (time.time() - start_time) / 60  # minutes
    print(f"Completed {sample_name} | Rows: {df.shape[0]} | Time: {elapsed:.2f} min")

counts_matrix = all_counts[0]
for df in all_counts[1:]:
    counts_matrix = counts_matrix.merge(df, on="Geneid", how="outer")

output_file = os.path.join(path, "GSE106305_counts_matrix_3011.csv")
counts_matrix.to_csv(output_file, index=False)

print("\All files processed!")
print("Merged matrix shape:", counts_matrix.shape)
print("Saved to:", output_file)

# Downstream Analysis — DESeq2 in R

This section covers:

* Importing count matrix

* Running DESeq2

* PCA / heatmap / volcano plot

* GO/KEGG enrichment

## Differential Expression Analysis (DESeq2) in R

Install required packages:

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

Run the R script - "DESeq2_tutorial_GSE106305_Updated.Rmd" @RStudio and generate all downstream statistical analysis

# End of Tutorial

This README provides a full, reproducible workflow for bulk RNA-seq—from raw FASTQ files to biological interpretation of differentially expressed genes.

## Acknowledgements

I would like to acknowledge **Smriti Arora** for conducting the workshop that guided the development of this RNA-seq analysis pipeline.
