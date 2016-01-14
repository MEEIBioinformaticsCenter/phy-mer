##Phy-Mer

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

####Citation:
Daniel Navarro-Gomez, Jeremy Leipzig, Lishuang Shen, Marie Lott, Alphons P.M. Stassen, Douglas C. Wallace, Janey L. Wiggs, Marni J. Falk, Mannis van Oven, and Xiaowu Gai
Phy-Mer: a novel alignment-free and reference-independent mitochondrial haplogroup classifier.

Bioinformatics (2015) 31 (8): 1310-1312. doi: 10.1093/bioinformatics/btu825  -> http://bioinformatics.oxfordjournals.org/content/31/8/1310



####1. Installation: 

#####1.1. Installing dependences:

Download and install the dependencies. Installation instructions for each dependency are provided on its project website with the URL provided.

- Python (2.7.3 and 3.0 tested -v3.0 thanks to @fmerinocasallo-) -> https://www.python.org/downloads/

- Pysam (0.7.4 tested) -> https://code.google.com/p/pysam/ 

- BioPython (1.58 tested) -> http://biopython.org/wiki/Download

#####1.2. Download Phy-Mer package, uncompress it and uncompress the Library:
Download last version of Phy-Mer and utilities following this link (https://github.com/danielnavarrogomez/phy-mer/archive/master.zip), uncompress it, open "phy-mer-master" folder and uncompress "PhyloTree_b16_k12.txt.gz".
	
	wget https://github.com/danielnavarrogomez/phy-mer/archive/master.zip
	unzip master.zip
	cd phy-mer-master
	gzip -d PhyloTree_b16_k12.txt.gz

####2. Usage

#####SYNOPSIS
	min_kmer_repeats=1
	./Phy-Mer.py [--verbose] [--print-ranking] [--def-snps=haplogroup_def_motifs.csv] [--min-DoC=10] Library.txt INPUT_1 [INPUT_2 ... INPUT_X]

#####DESCRIPTION
	
	Optinal arguments.
	  --help                     Show this help.
	  --verbose                  Print step by step process.
	  --print-ranking            Print top 5 results instead of the best match.
	  --def-snp=file.csv         Add Haplogroup defining snps to top matches
	                             based in file.csv (Build_16_-_rCRS-based_haplogroup_motifs.csv
	                             in resources folder).
	  --min-DoC=10               Only apply to BAM inputs. Minimal number of occurences of a K-mer
                                     to be consider.


	Mandatory arguments.
	  Library.txt                Provided with the package: PhyloTree_b16_k12.txt
	  INPUT_X                    One (required) or more input sequence files in FASTA, FASTQ, or BAM format.

####3. Phy-Mer Output
#####3.1. Normal output

	KF899911.fasta   I1a1a1  0.998459118106
 
**KF899911.fasta**:              Name of the input file
**I1a1a1**:                      Haplogroup prediction
**0.998459118106**:              Score (See Score section for more info)

#####3.2. Verbose output (*--verbose*)

	Opening Library and Checking k-mer size...
	READING WHOLE FILE IN MEMORY
	DONE
	K-mer=12
	Processing KF899911.fasta
	Opening fasta file and loading it in memory...
	Creating k-mer in memory (K=12)...
	16531 K-mers as input
	Comparing input with Library...
	Creating score table...
	KF899911.fasta   ['I1a1a1', 0.997856951513528, 0.9990612846989406, 0.9984591181062343]

#####3.3. Top 5 results output (*--print-ranking*)
	KF899911.fasta
	I1a1a1  0.998459118106
	I1a1a   0.997688614168
	I1a1a2  0.996884512927
	I1a1    0.996884512927
	I1a1e   0.996080411687

Top 5 scored results.

#####3.4. Haplogroup defining snps output (*--def-snp=file.csv*)
	KF899911.fasta   I1a1a1  0.998459118106  [['73', '199', '203', '204', '250', '263', '455.1T', '573.XC', '750', '1438', '1719', '2706', '3447', '3990', '4529T', '4769', '6734', '7028', '8251', '8616T', '8860', '9053', '9947', '10034', '10238', '10398', '10915', '11719', '12501', '12705', '13780', '14766', '15043', '15326', '15547', '15924', '16129', '16172', '16223', '16311', '16391']]


####4. Phy-Mer Scores

#####4.1. Standard score
Given a mitochondrial sequence as the input, Phy-Mer first decomposes its sequence into a set of all possible 12-mers, which are then compared against each of the 12-mer libraries of all haplogroups. As an example, an input sequence has x 12-mers that matches the 12-mers of haplogroup A, which has a total y 12-mers in its library. The sequence has on the other hand z 12-mers that match 12-mers of all haplogroup libraries. 

Then a score is derived as (x/y + x/z)/2.

#####4.2. Verbose Scores
Following previous section, we define verbose scores as:

[ 'SAMPLE.fasta', x/y , x/z , (x/y + x/z)/2 ]

####5. Utility Tools
#####5.1. Create a K-mer library:

######SYNOPSIS
	./build_Phy-Mer_Library.py MtDNA_Genome.fasta SNPS_HAPLOGROUPS.csv RESULT_Library.txt [min_k-mer_size]

######DESCRIPTION
**build_Phy-Mer_Library.py**: Python program that creates the k-mer library of all haplogroups.

**MtDNA_Genome.fasta**: Mitochondrial reference genome sequence in FASTA format, e.g. the revised Cambridge Reference Sequence (rCRS).

**SNPS_HAPLOGROUPS.csv**: The list of haplogroup-defining SNPs based on the mitochondrial reference genome sequence. An example file for PhyloTree Build 16 based in rCRS is provided in the package, named Build_16_-_rCRS-based_haplogroup_motifs.csv.

**RESULT_Library.txt**: The output k-mer library of all haplogroups defined in the Build_16_-_rCRS-based_haplogroup_motifs.csv file.

**min_k-mer_size**: OPTIONAL, minimal k-mer size considered for the library

######EXAMPLE
	Ex: ./build_Phy-Mer_Library.py MtGenome_sequence.fasta Build_16_-_rCRS-based_haplogroup_motifs.csv Custom_PhyloTree_b16.txt

#####5.2. Create a FASTA file from snp data:
######SYNOPSIS
	./convert_MtSNP_to_MtFasta.py REFERENCE.fasta SNP.csv RESULT.fasta
######DESCRIPTION
**convert_MtSNP_to_MtFasta.py**: Python script to create a FASTA formatted mtDNA sequence based on a list of SNPs

**REFERENCE.fasta**: mtDNA reference genome sequence upon which the SNPs are based

**SNP.csv**: SNP data in comm-delimited format.

**RESULT.fasta**: Output file to store the FASTA sequence


######EXAMPLE
	Ex: ./convert_MtSNP_to_MtFasta.py MtGenome_sequence.fasta example_snp.csv OUT.fasta



