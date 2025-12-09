import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path

# -----------------------------
# Configuration
# -----------------------------
# Define folder paths relative to the main directory
RAW_DATA_DIR = Path("Raw_Data")
FINAL_FASTA_DIR = Path("homologs")

SAMPLES = ["A", "B", "C", "D"]
parts = [1, 2, 3]

# Create output folder if it doesn't exist
if not FINAL_FASTA_DIR.exists():
    FINAL_FASTA_DIR.mkdir()

# -----------------------------
# Function to Clean FASTQ (Robust Mode)
# -----------------------------
def cleanup_fastq(filename):
    """
    Reads raw FASTQ, removes whitespace, and TRUNCATES sequence/quality
    to matching lengths to prevent errors.
    """
    temp_file = filename.with_suffix('.tmp.fq')
    
    with open(filename, 'r') as infile, open(temp_file, 'w') as outfile:
        while True:
            # Read 4 lines (one FASTQ record) manually
            header = infile.readline()
            if not header: break # End of file
            
            seq_line = infile.readline()
            plus_line = infile.readline()
            qual_line = infile.readline()
            
            # Clean whitespace from sequence and quality
            # We strip newlines first, then remove internal spaces
            clean_seq = seq_line.strip().replace(' ', '').replace('\t', '')
            clean_qual = qual_line.strip().replace(' ', '').replace('\t', '')
            
            # --- THE FIX: Force lengths to match ---
            if len(clean_seq) != len(clean_qual):
                # If they don't match, trim both to the shorter length
                min_len = min(len(clean_seq), len(clean_qual))
                clean_seq = clean_seq[:min_len]
                clean_qual = clean_qual[:min_len]
            
            # Write the fixed record to temp file
            outfile.write(header.strip() + '\n')
            outfile.write(clean_seq + '\n')
            outfile.write(plus_line.strip() + '\n')
            outfile.write(clean_qual + '\n')

    # Now BioPython can parse the corrected temp file safely
    for record in SeqIO.parse(temp_file, "fastq"):
        yield record
    
    # Clean up temp file
    if temp_file.exists():
        os.remove(temp_file)

print("Starting FASTQ cleaning...")

# -----------------------------
# Main Loop
# -----------------------------
for sample_id in SAMPLES:
    temp_records = []
    
    for part_num in parts:
        # Look for file inside Raw_Data folder
        raw_filename = RAW_DATA_DIR / f"sample{sample_id}_part{part_num}.FASTQ"
        
        if not raw_filename.exists():
            print(f"Warning: Missing {raw_filename}")
            continue
        
        try:
            for record in cleanup_fastq(raw_filename):
                # Normalize header id
                record.id = record.id.replace(f"part{part_num}", f"part_{part_num}")
                temp_records.append(record)
        except Exception as e:
            print(f"Error processing {raw_filename.name}: {e}")

    if temp_records:
        # Concatenate sequences
        final_seq_str = "".join([str(r.seq) for r in temp_records])
        final_record = SeqRecord(
            Seq(final_seq_str),
            id=f"sample_{sample_id}_full",
            description=f"Concatenated sequence Sample {sample_id}"
        )
        
        output_path = FINAL_FASTA_DIR / f"sample{sample_id}_final.fas"
        SeqIO.write(final_record, output_path, "fasta")
        print(f"Created: {output_path}")

print("Cleaning complete.")