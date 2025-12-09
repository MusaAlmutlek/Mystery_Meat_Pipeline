# Mystery_Meat_Pipeline
Analysis pipeline for the identification of mystery meat samples
Forensic DNA Barcoding Pipeline: Mystery Meat Analysis
Overview
This repository contains a bioinformatics pipeline developed to assist the Animal and Plant Health Agency (APHA) in a forensic investigation involving unidentified meat products intercepted at Heathrow Airport. The primary objective is to automate the processing of mitochondrial DNA (COI gene) sequencing data to identify the species of origin and determine if they belong to endangered or protected taxa.

Pipeline Workflow
The analysis is broken down into three automated stages:

Data Pre-processing: Cleaning raw FASTQ reads, quality filtering, and merging fragmented sequences.

Sequence Translation: Converting DNA sequences into amino acid sequences using the Vertebrate Mitochondrial genetic code.

Aggregation: Combining all processed samples and reference sequences into a master file ready for multiple sequence alignment (MSA) and phylogenetic tree construction.

Repository Contents & Usage
1. Scripts/clean_and_convert.sh
A Bash script designed to sanitize raw sequencing data.

Function: It standardizes file headers, converts FASTQ files to FASTA format, masks low-quality bases (Phred < 20), and concatenates split sequence parts into a single consensus sequence per sample.

Requirements:

Requires the seqtk toolkit (ensure it is installed and in your system PATH).

Input files must follow the naming convention: sampleX_partY (e.g., sampleA_part1).

Execute this script in the directory containing your raw data.

2. Scripts/DNA_translate.py
A Python script for robust biological translation.

Function: Iterates through DNA FASTA files, identifies the longest Open Reading Frame (ORF) in each sequence, and translates it into a protein sequence using Translation Table 2 (Vertebrate Mitochondrial).

Requirements:

Requires the Biopython library.

Input: Expects a directory named homologs containing the DNA FASTA files.

Output: Creates (or populates) a directory named protein_sequences.

3. Scripts/cat_and_clean.sh
A utility script for final dataset preparation.

Function: Aggregates individual protein FASTA files into a single master dataset (all_proteins.fasta). It also standardizes alignment formatting by replacing unknown amino acid residues ('X') with gap characters ('-').

Usage: Run this script after the translation step to generate the final input file for alignment tools like MUSCLE or ClustalW.

Dependencies & Installation
To replicate this analysis, ensure the following tools are installed on your system:

Python 3.x

Biopython Library:

Bash

pip install biopython
seqtk: Download and compile from lh3/seqtk on GitHub. Ensure the executable is accessible in your system's $PATH.

Reproducibility
All scripts are designed to be run sequentially. Please refer to the specific script descriptions above for input/output directory structures (Raw_Data, homologs, Results). Raw sequencing data is stored in the Raw_Data directory for transparency.

Author: Musa Husain G Almutlek (Student ID: 2627075) Date: December 2025
