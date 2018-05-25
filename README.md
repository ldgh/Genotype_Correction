# Genotype_Correction

This script was written to correct the genotyping errors of the GATK pipeline when genotyping an LDGH Haloplex Targeted Sequencing dataset. Many of the mistakes originated from low coverage regions called by GATK or from PCR duplicates.

# Requirements
The script is written in Python 3 and uses the pandas module to perform operations on the dataset.
1. A list of the alternative alleles in the loci of interest in a .tsv format which contains the chromosome, position and the alternative allele of the variant.
2. One or multiple .pileup files which also contain the information at which position in a read a base was called.

# Function
The script takes a list of variants that are meant to be genotyped and pileup files of samples in order to create vcf files containing the genotypes.

# Warning
The parameters of the script have to be adjusted to the dataset at hand, they are not universal, a different sequencing depth would require different parameters.
