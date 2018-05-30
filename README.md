# Genotype_Correction

This script was written in order to correct the genotyping errors of the GATK pipeline, when genotyping an LDGH Haloplex Targeted Sequencing dataset. Many of the errors of genotyping originated from low coverage regions called by GATK or from PCR duplicates. This tool with its parameters tailored to the dataset, can be utilized to genotype single nucleotide variants in a dataset that contains PCR duplicates (such as amplicon sequencing).

# Requirements
The script is written in Python 3 and uses the pandas module.
1. A list of the alternative alleles in the loci of interest in a .tsv format which contains the chromosome, position and the alternative allele of the variant. This list can be the list of variants called by the GATK pipeline.
2. One or multiple .pileup files which also contain the information at which position of the read a base was called (This latter information is neccessary for telling apart potential PCR duplicates and clearly independent reads).

# Function
The script takes a list of variants that are meant to be genotyped and pileup files of samples in order to create vcf-like .genotype files containing the genotypes. The coding of the genotypes follows the 0/1 format of .vcf files, however partial genotyping 1/. is also allowed if one allele can be called with certainity while the other cannot.

# Warning
The parameters of the script have to be adjusted to the dataset at hand, they are not universal! A different sequencing depth would require different parameters.
