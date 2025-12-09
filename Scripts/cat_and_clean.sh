#!/bin/bash

# Combine all protein sequences into one file
echo "Combining protein sequences..."
cat protein_sequences/*.fas > all_proteins.fasta

# Replace 'X' (unknown amino acid) with '-' (gap) for alignment
sed -i 's/X/-/g' all_proteins.fasta

echo "Done! Final file is: all_proteins.fasta"