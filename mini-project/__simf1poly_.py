import os 
import random #for random number gneration (will use for SNP positions)

"""
This script will load my parent1 and parent2 VCF files, and use them to simulate an F1 hybrid
by (re)combining alleles from the parental genomes.

The output file will be a VCF file containing the F1 genome. The parent VCF files should have 
REF, ALT, and GT fields. If they do not, you will want to go back and add those.

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
def parse_vcf(file_path):
    """
    Parse an existing VCF file into a list of variants.
    Each variant is a dictionary with CHROM, POS, REF, ALT, GENOTYPE.
    """
    variants = [] #empty list

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue  #skipping headers

            parts = line.strip().split('\t')

            chrom = parts[0]
            pos = parts[1]
            ref = parts[3]
            alt = parts[4]
            genotype = parts[9] if len(parts) > 9 else "./."

            variant = {
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "GT": genotype
            }

            variants.append(variant) #add variant to the list

    return variants #returns all variants

def sim_f1(parent1_variants, parent2_variants):
    """
    Combine two lists of parent variants to generate F1 hybrid variants.
    """
    #create a lookup table for each parent
    parent1_lookup = {}
    for var in parent1_variants:
        key = (var["CHROM"], var["POS"]) #look up a variant at a certain chromosome position
        parent1_lookup[key] = var

    parent2_lookup = {}
    for var in parent2_variants:
        key = (var["CHROM"], var["POS"])
        parent2_lookup[key] = var

    #gather unique variant positions (union of both parents)
    all_positions = set(parent1_lookup.keys()).union(set(parent2_lookup.keys()))

    f1_variants = [] #empty list

    for key in sorted(all_positions):
        chrom, pos = key

        #get parent variants
        var1 = parent1_lookup.get(key, None)
        var2 = parent2_lookup.get(key, None)

        ref = var1["REF"] if var1 else var2["REF"]
        alt1 = var1["ALT"] if var1 else ""
        alt2 = var2["ALT"] if var2 else ""

        # Build list of ALTs for F1
        alts = [] #empty list
        if alt1 != "":
            alts.extend(alt1.split(","))
        if alt2 != "":
            alts.extend(alt2.split(","))

        alts = list(sorted(set(alts)))  #remove any duplicate variants

        if not alts:
            continue  #skip if there are no variants

        #get genotypes from parents; assume 0/0
        gt1 = var1["GT"] if var1 else "0/0"
        gt2 = var2["GT"] if var2 else "0/0"

        #convert genotype strings like "0/1" to allele lists
        alleles1 = gt1.replace('|', '/').split('/')
        alleles2 = gt2.replace('|', '/').split('/')

        #random pick one allele from each parent
        allele_from_p1 = random.choice(alleles1)
        allele_from_p2 = random.choice(alleles2)

        #build F1 genotype
        f1_alleles = sorted([allele_from_p1, allele_from_p2])

        genotype = f"{f1_alleles[0]}/{f1_alleles[1]}"

        # Create the F1 variant record
        variant = {
            "CHROM": chrom,
            "POS": pos,
            "REF": ref,
            "ALT": ",".join(alts),
            "GT": genotype
        }

        f1_variants.append(variant)

    return f1_variants

def write_vcf(variants, output_file):
    """
    Write list of F1 variants to a VCF file
    """
    with open(output_file, 'w') as vcf:
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
    f1_variants = simulate_f1(parent1_variants, parent2_variants)

    # Step 4: Write the F1 hybrid VCF file
    write_vcf(f1_variants, f1vcf_output)

    print("[COMPLETE] F1 hybrid simulation finished!")









