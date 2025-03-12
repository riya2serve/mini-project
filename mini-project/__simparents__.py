from Bio import Entrez, SeqIO #BioPython for NCBI access and FASTA parsing
import os 
import random #for random number gneration (will use for SNP positions)
import pandas as pd #for dataframe handling (not necessary)

#================
#CONFIGURATION
#================

"""
This script will allow users to simulate two parental genomes from reference FASTA files.
It will introduce random SNPs to generate Parent 1 and Parent 2.

With these, the user can then simulate an F1 hybrid VCF. If SNPs differ, the F1 gets heterozygous (0/1);
if SNPs are identical the F1 gets homozygous (1/1 OR 0/0).
"""

genome_fasta = "spinach_genome/spinach_wg.fasta"

output_folder = "spinach_genome"
os.makedirs(output_folder, exist_ok = True) #checks that output folder exists

#================
#FUNCTION(s)
#================

def simulate_snps(fasta_file, snp_count, output_vcf):
    """
    Simulate random SNPs in a genome FASTA and save as a VCF.
    """
    sequences = list(SeqIO.parse(fasta_file, "fasta")) #can contain contigs or chromosomes
    tot_nome_length = sum(len(seq) for seq in sequences) #add up the contig or chromosome length so you can know full range for 'placing' SNPs

    print(f"[INFO] Loaded genome: {tot_nome_length:,} bp across {len(sequences)} chromosomes")

    # Generate random SNP positions
    snp_positions = random.sample(range(1, tot_nome_length), snp_count)
    snp_positions.sort()

    # Store SNPs per chromosome
    snps_on_chrom = {} #this is a dictionary for storing the SNPs on each chromsome
    genome_cursor = 0 #position on genome; initializing at zero

    for seq in sequences: #this for loop iterates through each chromosome
        chrom = seq.id #name of the chromosome
        chrom_length = len(seq)
        chrom_snps = [] #an empty list to store the SNPs found on a chromsome

        #will keep looping as long as these conditions are true:
        while snp_positions and genome_cursor < snp_positions[0] <= genome_cursor + chrom_length:
            pos_in_chrom = snp_positions.pop(0) - genome_cursor #converting genome wide SNPs into chrom-relative position
            ref_base = seq.seq[pos_in_chrom - 1].upper() #gets the reference nucleotide and converts it to uppercase

            if ref_base not in "ACGT":
                continue  # Skip ambiguous/unknown bases

            # Generate a random alternative allele
            alt_base = random.choice([b for b in "ACGT" if b != ref_base])
            chrom_snps.append((pos_in_chrom, ref_base, alt_base)) #add a tuple describing the SNP

        genome_cursor += chrom_length

        if chrom_snps:
            snps_on_chrom[chrom] = chrom_snps

    # Write VCF file
    with open(output_vcf, "w") as vcf:
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write(f"##source=Simulated_Parent\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPARENT\n")

        for chrom, snps in snps_on_chrom.items():
            for pos, ref, alt in snps:
                vcf.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\tGT\t1/1\n")

    print(f"[DONE] Simulated SNPs written to {output_vcf}")

# ================
# FUNC. EXECUTION
# ================
if __name__ == "__main__":

    # Simulate SNPs for Parent1
    simulate_snps(
        fasta_file=genome_fasta,
        snp_count=num_snps,
        output_vcf=os.path.join(output_folder, "parent1.vcf")
    )

    # Simulate SNPs for Parent2 (different SNP set)
    simulate_snps(
        fasta_file=genome_fasta,
        snp_count=num_snps,
        output_vcf=os.path.join(output_folder, "parent2.vcf")
    )

    print("[COMPLETE] Parental VCF simulation finished!")