#!/bin/bash

# Directory containing trimmed FASTQ files
FASTQ_DIR="trimmed"
GENOME_INDEX="grch38/genome"
LOGFILE="alignment_log.txt"

# Clear or create logfile
> "$LOGFILE"

# List of trimmed FASTQ files
FILES=(
    "LNCAP_Hypoxia_S1_trimmed.fastq.gz"
    "LNCAP_Hypoxia_S2_trimmed.fastq.gz"
    "LNCAP_Normoxia_S1_trimmed.fastq.gz"
    "LNCAP_Normoxia_S2_trimmed.fastq.gz"
    "PC3_Hypoxia_S1_trimmed.fastq.gz"
    "PC3_Hypoxia_S2_trimmed.fastq.gz"
    "PC3_Normoxia_S1_trimmed.fastq.gz"
    "PC3_Normoxia_S2_trimmed.fastq.gz"
)

for f in "${FILES[@]}"; do
    FILE_PATH="$FASTQ_DIR/$f"
    
    # Check if the file exists
    if [[ ! -f "$FILE_PATH" ]]; then
        echo "ERROR: File $FILE_PATH not found!" | tee -a "$LOGFILE"
        continue
    fi

    SAMPLE_NAME=$(basename "$f" _trimmed.fastq.gz)
    echo "Processing $SAMPLE_NAME at $(date)" | tee -a "$LOGFILE"

    # Run HISAT2 and pipe to samtools sort
    hisat2 -q -x "$GENOME_INDEX" -U "$FILE_PATH" 2>>"$LOGFILE" | \
    samtools sort -o "${SAMPLE_NAME}.bam"

    # Check if BAM was created
    if [[ -f "${SAMPLE_NAME}.bam" ]]; then
        samtools index "${SAMPLE_NAME}.bam"
        echo "Finished $SAMPLE_NAME at $(date)" | tee -a "$LOGFILE"
        echo "--------------------------------------" | tee -a "$LOGFILE"
    else
        echo "ERROR: Failed to create BAM for $SAMPLE_NAME" | tee -a "$LOGFILE"
    fi
done

echo "All alignments completed." | tee -a "$LOGFILE"

