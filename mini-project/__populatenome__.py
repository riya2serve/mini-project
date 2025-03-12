from Bio import Entrez, SeqIO
import os 
import random 
import pandas as pd 

#================
#CONFIGURATION
#================
#provide an email bc NCBI requires it?
Entrez.email = "rr3491@columbia.edu"

#================
#MY FUNCTIONS
#================

#Step 1. make an output directory
os.makedirs(output_folder, exist_ok = True)

output_FASTA = f"{output_folder}{accession_id}.fasta"


#let's create a function to fetch genome sequences from NCBI
def fetch_nomes(accession_id, output_folder = "spinach_genome/"):
    """
    This function will fetch genomesfrom NCBI using accession ids provided by the user. 
    Once located, the genome(s) will be saved as a FASTA file.

    Parameters:
    accession_id (string - str): NCBI accession number
    output_folder (str): folder to save the FASTA file(s) (default: current user director)
    """










#second genome is an 