##Phy-Mer

####How Phy-Mer works: 

#####One time step
Uncompress the Library:
	
	gzip -d PhyloTree_b16_k12.txt.gz

#####Usage
	./Phy-Mer.py [--verbose] [--print-ranking] [--def-snps=haplogroup_def_motifs.csv] Library.txt INPUT_1 [INPUT_2 ... INPUT_X]
	
	Optinal arguments.
	  --verbose                  Print step by step process.
	  --print-ranking            Print a ranking of best 5 results instead of the best match.
	  --def-snp=file.csv         Add Haplogroup defining snps based in file.csv (Build_16_-_rCRS-based_haplogroup_motifs.csv
	                             in resources folder) to the result.

	Mandatory arguments.
	  DataBase.txt               Provided with the package: PhyloTree_b16_k12.txt
	  INPUT_X                    Fasta, fastq, and bam files.


####How to create a Phy-Mer Library:

	./build_Phy-Mer_Library.py REFERENCE_FASTA_FILE.fasta SNPS_HAPLOGROUPS.csv RESULT_Library
	Ex: ./build_Phy-Mer_Library.py MtGenome_sequence.fasta Build_16_-_rCRS-based_haplogroup_motifs.csvv Custom_PhyloTree_b16.txt

####Tested with:
* Python 2.7.3
* Pysam 0.7.4
* BioPython 1.58
 
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


