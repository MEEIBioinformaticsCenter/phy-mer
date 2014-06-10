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
 


