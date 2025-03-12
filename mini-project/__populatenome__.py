from Bio import Entrez, SeqIO #BioPython for NCBI access and FASTA parsing
import os 
import random #for random number gneration (will use for SNP positions)
import pandas as pd #for dataframe handling (not necessary)

#================
#CONFIGURATION
#================
#provide an email bc NCBI requires it?
Entrez.email = "rr3491@columbia.edu"

#Step 1. make an output directory to store FASTA files
output_folder = "spinach_genome/"
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

#================
#MY FUNCTIONS
#================

#let's create a function to fetch genome sequences from NCBI
def fetch_fasta(accession_id, output_folder = "spinach_genome/"):
    """
    This function will fetch a FASTA file from NCBI using their accession IDs.
    Each FASTA will be saved locally (default: current user director)
    """
    print(f"[INFO] Fetching {accession_ids}. . .")

    try:
        #fetch sequences from NCBI
        handle = Entrez.efetch(db="nuccore", id=accession_ids, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()

        # Output file path
        output_path = os.path.join(output_folder, f"{accession_ids}.fasta")

        # Write FASTA(s) to output folder
        with open(output_path, "w") as fasta_file:
            fasta_file.write(fasta_data)

        print(f"[SUCCESS] Saved {accession_id} to {output_path}")

    except Exception as e:
        print(f"[ERROR] Failed to fetch {accession_id}: {e}")

# ==============
# FETCH CHROMS
# ==============
if __name__ == "__main__":
    for acc_id in accession_ids:
        fetch_fasta(acc_id, output_folder)

    print("[COMPLETE] Finished downloading all chromosomes.")







#second genome is an 