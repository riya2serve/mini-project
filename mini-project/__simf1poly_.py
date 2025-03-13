import os 
import random #for random number gneration (will use for SNP positions)
import argparse

"""
This script will load two parent VCF files and use them to simulate an F1 hybrid
by (re)combining alleles from the parental genomes. The output is an F1 hybrid VCF file that includes 
REF, ALT, and GT fields. 

Note: Parent VCF files should include REF, ALT, and GT fields. If they do not, you will want to go back 
and add them first.

"""

# ===============
# CONFIGURATION
# ===============
p1vcf_path = "spinach_genome/parent1.vcf"
p2vcf_path = "spinach_genome/parent2.vcf"

#output file and file path
f1vcf_output = "spinach_genome/f1hybrid.vcf"

# ============
# FUNCTION(s)
# ============

def get_alleles(gt):
    """
    Converts a genotype string (GT) into a list of allele indices.
    Arguments: 
        gt (string): genotype string
    Return:
        list: alleles expressed as integers
        "0/1" -> [0, 1]
        "./." -> [0] (cases where reference allele is missing)

    """
    if gt == "./.":
        return [0]  # assume reference allele, if the genotype is missing
    alleles = gt.replace('|', '/').split('/')
    return [int(allele) for allele in alleles]

def parse_vcf(file_path):
    """
    Parses an existing VCF file and extracts relevant information.
    Arguments:
        file_path (string): path to the parent VCF file
    Return:
        list: a list of dictionaries representing each variant.  
    """
    variants = [] #empty list to hold variant dictionaries; final will be list of dictionaries

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  #skipping headers

            #splits VCF lines by tab
            parts = line.strip().split('\t')

            #extracting specific VCF fields
            chrom = parts[0] #chromosome
            pos = parts[1] #position on chromosome
            ref = parts[3] #reference allele
            alt = parts[4] #alternate allele (s) -> polymorphic!
            genotype = parts[9] if len(parts) > 9 else "./."

            #formatting information that will be stored
            variant = {
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "GT": genotype
            }

            variants.append(variant) #adding to the list of variants

    return variants #returns all variants

def sim_f1(parent1_variants, parent2_variants):
    """
    Simulates F1 hybrid by combining parental variants
    Arguments:
        parent1_variants (list): variants from parent 1
        parent2_variants (list): variants from parent 2
    Return:
        list: simulated list of F1 variants
    """
    #creating dictionaries that can be looked up (CHROM,POS): variant
    parent1_lookup = {(v["CHROM"], v["POS"]): v for v in parent1_variants}
    parent2_lookup = {(v["CHROM"], v["POS"]): v for v in parent2_variants}

    #gathering all unique variant positions
    all_positions = set(parent1_lookup.keys()).union(set(parent2_lookup.keys()))

    f1_variants = [] #empty list that will hold all the f1 variants

    for chrom, pos in sorted(all_positions): #for-loop to iterate through variant positions
        #getting reference (REF) allele from both parents
        var1 = parent1_lookup.get((chrom, pos), None)
        var2 = parent2_lookup.get((chrom, pos), None)

        # Get REF allele from either parent
        ref = var1["REF"] if var1 else var2["REF"]

        #gathering all ALT alleles from parents
        alt_list = []
        if var1 and var1["ALT"] != ".":
            alt_list.extend(var1["ALT"].split(","))
        if var2 and var2["ALT"] != ".":
            alt_list.extend(var2["ALT"].split(","))

        #removing duplicates and sorting ALT alleles
        alts = sorted(set(alt_list))

        #if there are no ALT alleles, skip this site
        if not alts:
            continue

        #combines all ALT alleles into a comma-separated string!!!
        alt_field = ",".join(alts)

        #getting parent genotypes; if reference is missing use 0/0
        gt1 = var1["GT"] if var1 else "0/0"
        gt2 = var2["GT"] if var2 else "0/0"

        #get list of alleles from each parent
        alleles_p1 = get_alleles(gt1)
        alleles_p2 = get_alleles(gt2)

        # Randomly inherit one allele from each parent
        allele_from_p1 = random.choice(alleles_p1)
        allele_from_p2 = random.choice(alleles_p2)
        
        f1_alleles = sorted([allele_from_p1, allele_from_p2])
        
        #convert alleles into genotype strings for VCF format!!!
        genotype = f"{f1_alleles[0]}/{f1_alleles[1]}"

        if genotype == "0/0":
            continue

        # Create the variant record to be written
        variant = {
            "CHROM": chrom,
            "POS": pos,
            "REF": ref,
            "ALT": alt_field,
            "GT": genotype
        }

        f1_variants.append(variant)

    return f1_variants

def write_vcf(variants, output_file):
    """
    Write list of F1 variants to a VCF file
    """
    with open(f1vcf_output, 'w') as vcf:
        # Write VCF header
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##source=F1_hybrid_simple_simulation\n")
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF1_Hybrid\n")

        # Write each variant
        for var in variants:
            vcf.write(f"{var['CHROM']}\t{var['POS']}\t.\t{var['REF']}\t{var['ALT']}\t.\tPASS\t.\tGT\t{var['GT']}\n")

    print(f"[DONE] F1 hybrid VCF saved to: {f1vcf_output}")

# ==========================================
# MAIN PROGRAM
# ==========================================
if __name__ == "__main__":

    # Step 1: Check files exist
    if not os.path.exists(p1vcf_path) or not os.path.exists(p2vcf_path):
        print("[ERROR] One or both parent VCF files are missing!")
        exit(1)

    # Step 2: Parse parent VCF files
    print("[INFO] Reading Parent1 VCF...")
    parent1_variants = parse_vcf(p1vcf_path)

    print("[INFO] Reading Parent2 VCF...")
    parent2_variants = parse_vcf(p2vcf_path)

    # Step 3: Simulate the F1 hybrid
    print("[INFO] Simulating F1 hybrid...")
    f1_variants = sim_f1(parent1_variants, parent2_variants)

    # Step 4: Write the F1 hybrid VCF file
    write_vcf(f1_variants, f1vcf_output)

    print("[COMPLETE] F1 hybrid simulation finished!")









