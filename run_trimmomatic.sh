#!/bin/bash

FASTQ_DIR="fastq"
TRIM_DIR="trimmed"
TRIMMOMATIC_JAR="Trimmomatic-0.39/trimmomatic-0.39.jar"

# Create trimmed output folder
mkdir -p $TRIM_DIR

# List of sample names
SAMPLES=(
"LNCAP_Hypoxia_S1"
"LNCAP_Hypoxia_S2"
"LNCAP_Normoxia_S1"
"LNCAP_Normoxia_S2"
"PC3_Hypoxia_S1"
"PC3_Hypoxia_S2"
"PC3_Normoxia_S1"
"PC3_Normoxia_S2"
)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Running Trimmomatic for: $SAMPLE"

    java -jar $TRIMMOMATIC_JAR SE \
        -threads 4 \
        $FASTQ_DIR/${SAMPLE}.fastq.gz \
        $TRIM_DIR/${SAMPLE}_trimmed.fastq.gz \
        TRAILING:10 -phred33

    echo "Done: $SAMPLE"
done

echo "All samples trimmed successfully."
