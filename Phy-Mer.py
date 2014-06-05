#!/usr/bin/env python

import sys
from Bio import SeqIO
from operator import itemgetter
import ast
import os
import pysam
from Bio import SeqIO, Seq, SeqRecord


#### FUNCTIONS

# Converting a fasta sequence (string) to a dictionary (hash table) of k-mers
def convert_to_k_mer_hash(array_seq,min_repeats=1):
	global K_MER_SIZE
	return_k_mer={}
	for seq in array_seq:
		count_base=0
		while count_base<= len(seq)-K_MER_SIZE:
			coun_kmer=0
			kmer=''
			while coun_kmer<K_MER_SIZE:
				kmer+=seq[count_base+coun_kmer]
				coun_kmer+=1
			try:	
				return_k_mer[kmer]+=1
			except KeyError:
				return_k_mer[kmer]=1
			count_base+=1
	return_k_mer_2={}
	for i in return_k_mer.keys():
		if return_k_mer[i]>=min_repeats:
			return_k_mer_2[i]=True
	return return_k_mer_2

# Reading variables from file hash (by lines)'s line
def read_var_from_file_line(file_name,line):
	global REF_INDEX_ARRAY
	to_return=[]
	numbers_array=file_name[line-1].split(', ')
	for i in numbers_array:
		to_return.append(REF_INDEX_ARRAY[int(i)])
	return to_return
	

# Generator to convert BAM files into Biopython SeqRecords.
def bam_to_rec(in_file):
	MT_names=['M','MT','chrM','chrMT','chrMt','Mt','ChrM','ChrMT','ChrMt','CHRM','CHRMT','CHRMt']
	bam_file = pysam.Samfile(in_file, "rb")
	while not('region' in locals()):
		try:
			region=bam_file.fetch(MT_names[0])
		except ValueError:
			MT_names.pop(0)
			if len(MT_names)==0:
				print "ERROR: MtDNA not found in your bam file: ['M','MT','chrM','chrMT','chrMt','Mt','ChrM','ChrMT','ChrMt','CHRM','CHRMT','CHRMt']"
				exit(1)
			pass
	for read in region:
		seq = Seq.Seq(read.seq)
		if read.is_reverse:
			seq = seq.reverse_complement()
		rec = SeqRecord.SeqRecord(seq, read.qname, "", "")
		yield rec

#### MAIN
def main():
	if len(sys.argv)!=3:
		print "ERROR: Usage: "+str(sys.argv[0])+" DataBase.txt INPUT_FASTA.fasta"
		exit(1)

	min_kmer_repeats=1
	FASTA_FILE=sys.argv[2]
	DB_FILE=sys.argv[1]
	
	print "Openning DB and Checking k-mer size..."	
	print "READING WHOLE FILE IN MEMORY"
	with open(DB_FILE) as f:
		TEST_content = f.readlines()
	print "DONE"
	
	
	REF_INDEX_ARRAY=ast.literal_eval(TEST_content[2])
		
	k_mer_hash=ast.literal_eval(TEST_content[0])
	K_MER_SIZE=len(k_mer_hash.keys()[1])
	print "K-mer="+str(K_MER_SIZE)
	
	total_hits=0
	result_haplogroup_ranking={}
		
	print "Openning fasta file and loading it in memory..."
	try:
		handle = open(FASTA_FILE, "rU")
	except IOError:
		print FASTA_FILE+" doesn't exist..."
		exit(1)
	array_seq=[]
	try:
		for record in SeqIO.parse(handle, "fasta") :
			array_seq.append(list(record.seq))
		handle.close()
		if len(array_seq)==0:
			print "Warning, no FASTA Reads, trying with FASTQ..."
			try:
				handle = open(FASTA_FILE, "rU")
			except IOError:
				print FASTA_FILE+" doesn't exist..."
				exit(1)	
			array_seq=[]
			for record in SeqIO.parse(handle, "fastq") :
				array_seq.append(list(record.seq))
			handle.close()
	except IndexError:
		min_kmer_repeats=10
		array_seq=[]
		print "Seems that is not a FASTA/FASTQ file... Trying with BAM..."
		for i in bam_to_rec(FASTA_FILE):
			array_seq.append(str(i.seq))
		pass	
		
	if len(array_seq)==0:
		print "ERROR: empty input"
		exit(1)

	print "Creating k-mer in memory (K="+str(K_MER_SIZE)+")..."
	input_fasta_k_mer=convert_to_k_mer_hash(array_seq,min_kmer_repeats)
	print str(len(input_fasta_k_mer))+" K-mers as input"
	print "Comparing input with DB..."
	total_hits=0
	
	for input_k_mer in input_fasta_k_mer.keys():
		try:
			db_line=k_mer_hash[input_k_mer]
			aux_array_ref=read_var_from_file_line(TEST_content,db_line)
			for aux_ref in aux_array_ref:
				try:
					result_haplogroup_ranking[aux_ref]+=1
				except KeyError:
					result_haplogroup_ranking[aux_ref]=1
					pass
			total_hits+=1
		except KeyError:
			pass
	print "Creating score table..."
        max_haplogroup_score_table=ast.literal_eval(TEST_content[1])

	ranking_table=[]
	for haplogroup in result_haplogroup_ranking.keys():
		ranking_table.append([haplogroup,float(result_haplogroup_ranking[haplogroup])/float(max_haplogroup_score_table[haplogroup]),float(result_haplogroup_ranking[haplogroup])/float(total_hits),((float(result_haplogroup_ranking[haplogroup])/float(max_haplogroup_score_table[haplogroup]))+(float(result_haplogroup_ranking[haplogroup])/float(total_hits)))/2])
	ranking_table.sort(key=itemgetter(3),reverse=True)

	i=0
	score=0.00
	result=['',0.00,0.00]
	try:
		while i<10:  ## we are printing only top scores
			if score>ranking_table[i][3]:
				break
			if result[0]=='':
				result=ranking_table[i]
			else:
				result[0]+=", "+ranking_table[i][0]
				
			score=ranking_table[i][3]
			i+=1
	except IndexError:
		print "ERROR: no result, check input file"
		exit(1)
	print result

##############
if __name__ == "__main__":
	main()


