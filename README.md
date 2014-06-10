#### How Phy-Mer works 

Uncompress the DB:
	
	gzip -d PhyloTree_b16_k12.txt.gz

Enjoy Phy-Mer:
	
	./Phy-Mer.py PhyloTree_b16_k12.txt INPUT*
	*INPUT -> fasta, fastq and bam compatible.


#### How to create a Phy-Mer DB 

	./build_Phy-Mer_DB.py REFERENCE_FASTA_FILE.fasta SNPS_HAPLOGROUPS.csv RESULT_DB


