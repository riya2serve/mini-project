from Bio import Entrez, SeqIO

#provide an email bc NCBI requires it
Entrez.email = "rr3491@columbia.edu"

#let's create a function to fetch genome sequences from NCBI
def fetch_nomes(accession_ids, output_folder = ""):
    """
    This function will fetch genome sequences from NCBI using the accession numbers provided by the user.
    Once located, the acessions will be saved as a FASTA file.

    accession_ids (parameter): user to provide NCBI accession number
    output_folder (parameter): folder to save the FASTA files (default: current user director)
    """

    for accession_id in accession_ids:
        output_fasta = f"{output_folder}{accession_id}.fasta"

        try:
            #fetch sequence(s) in FASTA format
            genome = Entrez.efetch(db = "nucleotide", id = accession_id, rettype = "fasta", retmode = "text")

            #save to output file
            with open(output_fasta, "w") as file:
                file.write(handle.read())

            handle.close()
            print(f"Genome {accession_id} saved as {output_fasta}")

        except Exception as e:
            print(f"Error in fetching genome {accession_id}: {e}")

fetch_nomes(["GCA_900095335.1", ""]) 
#first genome is an S. latifolia (male)
#second genome is an 