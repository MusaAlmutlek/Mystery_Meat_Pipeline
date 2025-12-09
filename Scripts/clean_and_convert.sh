#!/bin/bash

# ----------------------------------------------------
# MYSTERY MEAT ANALYSIS PIPELINE - Data Cleaning
# ----------------------------------------------------
# This script processes raw FASTQ files for Samples A-D.
# It performs the following steps:
# 1. Corrects formatting errors in headers.
# 2. Removes duplicate entries (Sample B).
# 3. Converts FASTQ to FASTA using seqtk (masking low quality bases).
# 4. Concatenates split parts into a single file per sample.
# ----------------------------------------------------

# Safety settings: Exit on error
set -e

# Check if seqtk is installed
if ! command -v seqtk &> /dev/null; then
    echo "Error: seqtk is not installed or not in your PATH."
    exit 1
fi

echo "Starting cleaning pipeline..."

# ----------------------------------------------------
# MAIN PROCESSING LOOP
# ----------------------------------------------------
# We loop through all 4 samples (A, B, C, D)
for sample in A B C D; do
    
    echo "Processing Sample $sample..."
    
    # Create an empty variable to store the names of the converted files
    converted_parts=""

    # Loop through the 3 parts for this sample
    for part in 1 2 3; do
        
        # Define filenames
        raw_file="sample${sample}_part${part}.FASTQ"
        clean_fasta="sample${sample}_part${part}.fasta"
        
        # Check if file exists
        if [[ ! -f "$raw_file" ]]; then
            echo "Warning: $raw_file not found. Skipping."
            continue
        fi

        # --- STEP 1: FIX HEADERS ---
        # The raw data has inconsistent headers.
        # We use sed to insert the underscore if it's missing.
        sed "s/sample${sample}_part${part}/sample${sample}_part_${part}/" "$raw_file" > temp_header.fq
        mv temp_header.fq "$raw_file"

        # --- STEP 2: REMOVE DUPLICATES (Specific to Sample B Part 1) ---
        if [[ "$sample" == "B" && "$part" == "1" ]]; then
            echo "  -> Removing duplicate header in Sample B Part 1..."
            awk '!(($0 == prev) && /^@/) {print} {prev = $0}' "$raw_file" > temp_dedup.fq
            mv temp_dedup.fq "$raw_file"
        fi

        # --- STEP 3: CONVERT TO FASTA ---
        # Use seqtk to convert to FASTA and mask low quality bases (q<20)
        echo "  -> Converting Part $part..."
        seqtk seq -a -q 20 -n N "$raw_file" > "$clean_fasta"

        # Add this new fasta file to our list
        converted_parts="$converted_parts $clean_fasta"
    done

    # --- STEP 4: CONCATENATE & FINALIZE ---
    # Merge all parts into one file
    cat $converted_parts > "sample${sample}_combined.fasta"

    # Create the final clean file with a single clean header
    echo ">sample_${sample}_full" > "sample${sample}_final.fas"
    grep -v "^>" "sample${sample}_combined.fasta" | tr -d '\n' >> "sample${sample}_final.fas"
    
    # Add a new line at the end
    echo "" >> "sample${sample}_final.fas"

    echo "Sample $sample complete. Final file: sample${sample}_final.fas"

done

# ----------------------------------------------------
# CLEANUP
# ----------------------------------------------------
echo "Cleaning up intermediate files..."
rm -f *_part*.fasta *_combined.fasta

echo "Pipeline completed successfully."
