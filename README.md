##Phy-Mer

####1. Installation: 

#####1.1. Installing dependences:

Python (2.7.3 tested) -> https://www.python.org/downloads/

Pysam (0.7.4 tested) -> https://code.google.com/p/pysam/ 

BioPython (1.58 tested) -> http://biopython.org/wiki/Download

#####1.2. Uncompress the Library
This is a one time step.
	
	gzip -d PhyloTree_b16_k12.txt.gz

####2. Usage

#####SYNOPSIS
	./Phy-Mer.py [--verbose] [--print-ranking] [--def-snps=haplogroup_def_motifs.csv] Library.txt INPUT_1 [INPUT_2 ... INPUT_X]

#####DESCRIPTION
	
	Optinal arguments.
	  --verbose                  Print step by step process.
	  --print-ranking            Print a ranking of best 5 results instead of the best match.
	  --def-snp=file.csv         Add Haplogroup defining snps based in file.csv (Build_16_-_rCRS-based_haplogroup_motifs.csv
	                             in resources folder) to the result.

	Mandatory arguments.
	  DataBase.txt               Provided with the package: PhyloTree_b16_k12.txt
	  INPUT_X                    Fasta, fastq, and bam files.

####3. Utility Tools
#####3.1. Create a K-mer Library:

######SYNOPSIS
	./build_Phy-Mer_Library.py MtDNA_Genome.fasta SNPS_HAPLOGROUPS.csv RESULT_Library.txt

######DESCRIPTION
build_Phy-Mer_Library.py: Python program that creates the k-mer library of all haplogroups.

MtDNA_Genome.fasta: Mitochondrial reference genome sequence in FASTA format, e.g. the revised Cambridge Reference Sequence (rCRS).

SNPS_HAPLOGROUPS.csv: The list of haplogroup-defining SNPs based on the mitochondrial reference genome secuence. An example file for PhyloTree Build 16 based in rCRS is provided in the package, named Build_16_-_rCRS-based_haplogroup_motifs.csv.

RESULT_Library.txt: The output k-mer library of all haplogroups defined in the Build_16_-_rCRS-based_haplogroup_motifs.csv file.

######EXAMPLE
	Ex: ./build_Phy-Mer_Library.py MtGenome_sequence.fasta Build_16_-_rCRS-based_haplogroup_motifs.csv Custom_PhyloTree_b16.txt

####3.2. Create a FASTA file from snp data:
######SYNOPSIS
	./convert_MtSNP_to_MtFasta.py REFERENCE.fasta SNP.csv RESULT.fasta
######DESCRIPTION
convert_MtSNP_to_MtFasta.py: Python script to create a FASTA formatted mtDNA sequence based on a list of SNPs

REFERENCE.fasta: mtDNA reference genome sequence upon which the SNPs are based

SNP.csv: SNP data in comm-delimited format.

RESULT.fasta: Output file to store the FASTA sequence


######EXAMPLE
	Ex: ./convert_MtSNP_to_MtFasta.py MtGenome_sequence.fasta example_snp.csv OUT.fasta

####Licensing:
    Phy-Mer
    Copyright (C) 2014  Daniel Navarro-Gomez (Daniel_navarro@meei.harvard.edu)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.


