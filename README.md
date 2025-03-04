
 my first mini-project

# VCF Simulation

This mini-project is designed to simulate a Variant Call Format (VCF) file that could be manipulated for genetic analysis.

The mock VCF file will contain a chromosome map for a single individual. The program will allow users to:

- Input two genome files and generate simulated genetic crosses
- Apply segregation distortion patterns to create mock F1 chromosomes
- Simulate chromosomes that are 1Mb in length, each with 1000 SNPs
 
The program will accept two parental genome files (**FASTA** or **VCF** format) and user-defined parameters for crossing and segregation patterns.  

The output will consist of:
- A simulated VCF file representing the recombined F1 chromosomes
- A summary of the genetic variations and **segregation distortion patterns** applied.

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
conda install [list dependencies here...] -c conda-forge
```
Alternatively, you can install via 'pip':
```
pip install -e .
```
4. Once installed, execute the program using
```
python __main__.py
```

### Using the CLI Tool

The program can be run from the terminal command line:
```bash
python __main__.py --parent1 parent1.vcf --parent2 parent2.vcf -- output simulated_F1.vcf --snp-count 1000 --chrom-length 1000000
```
```bash
usage: python __main__.py [-h] --parent1 PARENT1 --parent2 PARENT2 --output OUTPUT [--snp-count SNP_COUNT] [--chrom-length CHROM_LENGTH] [--distortion]

optional arguments:
 -h, --help			show this help message and exit 
 --parent1 PARENT1 		path to first parental genome file (FASTA/VCF)
 --parent2 PARENT2 		path to second parental genome file (FASTA/VCF)
 --output OUTPUT 		path to the simulated F1 chromosome VCF file
 --snp-count SNP_COUNT  number of SNPs to simulate per chromosome (default: 1000)
 --chrom-length CHROM_LENGTH  chromosome length in base pairs (default: 1Mb)
 --distortion 			apply segregation distortion 
```

These arguments will allow you to:
- Simulate recombination using two parental genomes
- Specify where the new F1 VCF file should be saved 
- Alter the simulation, specifically the # of SNPs and length of each chromosome
- Toggle with segregation distortion should you choose to implement it as a feature

### Sample Output

After running the program, the output of the VCF file should look something like this:
```simulated_F1.vcf
##fileformat = VCFv4.2
##source = SimulatedVCF
##contig =<ID=1, length = 1000000>

#CHROM 		POS 	ID 		REF 	ALT 	QUAL 	FILTER 	INFO
1 		1023    .       	A       G       .       .       DP = 13 #depth of coverage
1 		2071    .      		C       T       .       .       DP = 28
1 		3200    .       	C       T       .       .       DP = 22
```
