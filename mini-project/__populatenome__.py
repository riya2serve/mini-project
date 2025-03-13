from Bio import Entrez, SeqIO #BioPython for NCBI access and FASTA parsing
import os 
import random #for random number gneration (will use for SNP positions)
import pandas as pd #for dataframe handling (not necessary)

# ================
# CONFIGURATION
# ================
#requirement of NCBI entrez tools 
Entrez.email = "rr3491@columbia.edu"

#Step 1. define an output directory to store the downloaded FASTA files
output_folder = "spinach_genome"
#creates the directory if it doesn't exist already
os.makedirs(output_folder, exist_ok = True)

#Step 2. list the RefSeq chromosome accession IDs (on NCBI)
accession_ids = [
    "NC_079487.1",  # Chromosome 1
    "NC_079488.1",  # Chromosome 2
    "NC_079489.1",  # Chromosome 3
    "NC_079490.1",  # Chromosome 4
    "NC_079491.1",  # Chromosome 5
    "NC_079492.1"   # Chromosome 6
    #can add the mitochondrial (MT) and plastid (Pltd) IDs as needed
]

# ================
# MY FUNCTIONS
# ================

def fetch_fasta(accession_ids, output_folder):
    """
    This function will fetch a FASTA file from NCBI using their accession IDs.
    Each FASTA will be saved locally (default: current user director)
    
    parameters:
        accession_ids (string): NCBI accession number of a chromosome or nucleotide sequence
        output_folder (string): directory where downloaded FASTA files will be stored
    """
    print(f"[INFO] Fetching {accession_ids}. . .")

    try:
        #use Entrez.efectch to reqtireve sequences, from nuccore db on NCBI
        #fetch sequences from NCBI
        handle = Entrez.efetch(
            db="nuccore", 
            id=accession_ids, #accession IDs of chromosomes
            rettype="fasta", 
            retmode="text" #provides a plain text FASTA file
        ) 

        fasta_data = handle.read() #read FASTA data
        handle.close()

        #output file path
        output_path = os.path.join(output_folder, f"{accession_ids}.fasta")

        #open file in 'write mode' and save FASTA data
        with open(output_path, "w") as fasta_file:
            fasta_file.write(fasta_data)

        print(f"[SUCCESS] Saved {accession_ids} to {output_path}")

#if (above) doesn't work, then  print an error message and exception
    except Exception as e:
        print(f"[ERROR] Failed to fetch {accession_ids}: {e}")

# ==============
# FETCH CHROMS
# ==============
if __name__ == "__main__":
    for acc_id in accession_ids:
        fetch_fasta(acc_id, output_folder)

    print("[COMPLETE] Finished downloading all chromosomes.")

"""
==================================
EXPECTED OUTPUT (in terminal)
==================================

[INFO] Fetching NC_079487.1 from NCBI...
[SUCCESS] Saved NC_079487.1 to spinach_genome/NC_079487.1.fasta
"""










