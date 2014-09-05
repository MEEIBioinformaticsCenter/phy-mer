##Phy-Mer

####How Phy-Mer works: 

Uncompress the DB:
	
	gzip -d PhyloTree_b16_k12.txt.gz

Enjoy Phy-Mer:
	
	./Phy-Mer.py PhyloTree_b16_k12.txt INPUT*
	*INPUT -> fasta, fastq and bam compatible.


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


