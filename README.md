##Phy-Mer

####How Phy-Mer works: 

#####One time step
Uncompress the DB:
	
	gzip -d PhyloTree_b16_k12.txt.gz

#####Usage
	../Phy-Mer_src/Phy-Mer.py [--verbose] [--print-ranking] [--def-snps=haplogroup_def_motifs.csv] DataBase.txt INPUT_1 [INPUT_2 ... INPUT_X]
	Novel mitochondrial genome haplogroup defining algorithm using a k-mer approach.
	
	Optinal arguments.
	  --verbose                  Print step by step process.
	  --print-ranking            Print a ranking of best 5 results instead of the best match.
	  --def-snp=file.csv         Add Haplogroup defining snps based in file.csv (Build_16_-_rCRS-based_haplogroup_motifs.csv
	                             in resources folder) to the result.


####How to create a Phy-Mer DB:

	./build_Phy-Mer_DB.py REFERENCE_FASTA_FILE.fasta SNPS_HAPLOGROUPS.csv RESULT_DB
	Ex: ./build_Phy-Mer_DB.py resources/MtGenome_sequence.fasta resources/Build\ 16\ -\ rCRS-based\ haplogroup\ motifs.csv Custom_PhyloTree_b16.txt

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


