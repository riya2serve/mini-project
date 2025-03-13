
Spinach F1 Hybrid Simulation Project

# VCF Simulation

This mini-project simulates a Variant Call Format (VCF) file representing an F1 hybrid genome, which can be manipulated for genetic analysis and variant studies.

The mock VCF file will contain a chromosome map for a single individual. The program will allow users to:

- Input two genome files and generate simulated genetic crosses
- Apply segregation distortion patterns to create mock F1 chromosomes
- Simulate chromosomes that are 1Mb in length, each with 1000 SNPs
 
The program accepts two parental variant files (**FASTA** or **VCF** format) and user-defined parameters for simulating recombination and SNP inheritance. 

The output will consist of:
- A simulated VCF file representing the F1 hybrid chromosomes (combined variants from both parents)
- [planned] **segregation distortion patterns** and summaries of inheritance

### Installation

To install and run the program locally, follow these steps:

1. Clone the repository to your local computer and populate your new repo
```sh
cd myproject-name #open repo directory
git clone [https://github.com/YOUR-USERNAME/YOUR-REPOSITORY] #clone my github repo
```
After you do that you'' want to create a README file and some text to it to describe [...]

2. Install the required dependencies. If using **Conda**, you can install the required dependencies with:
```
conda install python=3.10 biopython -c conda-forge
```
Alternatively, you can install via 'pip':
```
pip install biopython
```
4. Once installed, execute the program using

### Using the CLI Tool

The program can be run from the terminal command line:
```bash
python __simparents__.py --parent1 spinach_genome/parent1.vcf --parent2 spinach_genome/parent2.vcf --snp-count 1000

python __simf1poly__.py
```
```bash
usage: python __simparents__.py
       --parent1 PARENT1
       --parent2 PARENT2
       [--snp-count 1000]
       [--chrom-length 1000000]
       [--distortion]

usage: python __simf1poly__.py
       --parent1 PARENT1
       --parent2 PARENT2
       --output OUTPUT

optional arguments:
 -h, --help			show this help message and exit
 --parent1 PARENT1 		path to first parental file (FASTA/VCF)
 --parent2 PARENT2 		path to second parental file (FASTA/VCF)
 --output OUTPUT 		path to the simulated F1 chromosome VCF file
 --snp-count SNP_COUNT  number of SNPs to simulate per chromosome (default: 1000)
 --chrom-length CHROM_LENGTH  chromosome length in base pairs (default: 1Mb)
 --distortion 			apply segregation distortion 
```

These arguments will allow you to:
- Simulate recombination using two parental variant files
- Specify where the new F1 VCF file should be saved 
- Alter the simulation, specifically the # of SNPs and length of each chromosome
- Toggle with egregation distortion [in progress]

### Sample Output

After running the program, the output of the VCF file should look something like this:
```f1hybrid.vcf
##fileformat = VCFv4.2
##source = _simparents_.py
##contig =<ID=1, length = 1000000>

#CHROM  POS   ID   REF  ALT  QUAL  FILTER  INFO  FORMAT  F1_Hybrid
1       1093  .    A    G    .     PASS    .     GT      0/1
1       2011  .    C    T    .     PASS    .     GT      1/1
```
