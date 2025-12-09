import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path 

# -----------------------------
# Configuration
# -----------------------------
# Use Path objects for reliable directory access on Windows
input_dir = Path("homologs")
output_dir = Path("protein_sequences")

# Create the output directory if it doesn't exist
if not output_dir.exists():
    output_dir.mkdir()

print(f"Reading files from {input_dir}...")

# -----------------------------
# Process each FASTA file
# -----------------------------

# We need to find BOTH .fas (your samples) and .fasta (your references) files
input_files = list(input_dir.glob("*.fas")) + list(input_dir.glob("*.fasta"))

for fasta_file in input_files: 
    
    # Define output path (keeps the same name but adds _prot)
    output_file = output_dir / f"{fasta_file.stem}_prot.fas"
    
    translated_records = []

    # Read sequences 
    # Use 'try' block to skip bad files without crashing
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            dna_seq = record.seq
            
            # Track the longest ORF found across all 3 frames
            best_orf = ""
            best_frame = 0
            
            # Check all three forward reading frames (0, 1, and 2)
            for frame in range(3):
                # Translate using Vertebrate Mitochondrial code (Table 2)
                # to_stop=False means we translate everything first
                protein = dna_seq[frame:].translate(table=2, to_stop=False)

                # Split at stop codons ("*") to find valid open reading frames
                fragments = str(protein).split("*")
                longest_in_frame = max(fragments, key=len)

                # If this is the longest ORF so far, save it
                if len(longest_in_frame) > len(best_orf):
                    best_orf = longest_in_frame
                    best_frame = frame

            # Create a new SeqRecord for the protein
            protein_seq = Seq(best_orf)
            
            new_record = SeqRecord(
                protein_seq,
                id=record.id,
                description=record.description,
                annotations={"frame": best_frame}
            )
            
            translated_records.append(new_record)

        # Write the translated proteins to the output folder
        if translated_records:
            SeqIO.write(translated_records, output_file, "fasta")
            print(f" -> Translated {fasta_file.name}")
        else:
            print(f" Warning: No sequences found in {fasta_file.name}")

    except Exception as e:
        print(f" Error processing {fasta_file.name}: {e}")

print("Translation complete.")