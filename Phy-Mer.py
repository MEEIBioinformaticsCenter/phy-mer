#!/usr/bin/env python

import sys
#import getopt
#import commands
#import string
from Bio import SeqIO
from operator import itemgetter
#from operator import attrgetter
import ast


import os
import pysam
from Bio import SeqIO, Seq, SeqRecord




#### FUNCTIONS

def convert_to_k_mer_hash(array_seq,min_repeats=1):
	global K_MER_SIZE
	return_k_mer={}
	#print "Starting with "+haplogroup_name
	#array_seq=reconstruct_array_seq(haplogroup_fasta[haplogroup_name])
	#print array_seq[0:15]
	#array_seq
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

def read_var_from_file_line(file_name,line):
	global REF_INDEX_ARRAY
	#return ast.literal_eval(linecache.getline(DB_FILE,db_line))
	#return ast.literal_eval(file_name[line-1])
	#return eval(file_name[line-1])
	to_return=[]
	numbers_array=file_name[line-1].split(', ')
	for i in numbers_array:
		to_return.append(REF_INDEX_ARRAY[int(i)])
	return to_return
	

def bam_to_rec(in_file):
	#"""Generator to convert BAM files into Biopython SeqRecords."""
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
	#print region
	#for read in bam_file:
	for read in region:
		
		#print read
		seq = Seq.Seq(read.seq)
		if read.is_reverse:
			seq = seq.reverse_complement()
		rec = SeqRecord.SeqRecord(seq, read.qname, "", "")
		yield rec



#### MAIN

if len(sys.argv)==3:
	min_kmer_repeats=1
	try:
		FASTA_FILE=sys.argv[2]
	except IndexError:
		FASTA_FILE="--interactive--"
	DB_FILE=sys.argv[1]
	
	print "Openning DB and Checking k-mer size..."	
	print "READING WHOLE FILE IN MEMORY"
	with open(DB_FILE) as f:
		TEST_content = f.readlines()
	print "DONE"
	
	
	REF_INDEX_ARRAY=ast.literal_eval(TEST_content[2])
		
	#k_mer_hash=ast.literal_eval(linecache.getline(DB_FILE,1))
	k_mer_hash=ast.literal_eval(TEST_content[0])
	K_MER_SIZE=len(k_mer_hash.keys()[1])
	print "K-mer="+str(K_MER_SIZE)
	#k_mer_hash=convert_all_hash_to_number_hash(k_mer_hash)
	
	
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
		#from Bio import AlignIO
		#from Bio.Align import AlignInfo
		#from Bio.Alphabet import IUPAC, Gapped
		#alphabet = Gapped(IUPAC.ambiguous_dna)
		print "Seems that is not a FASTA/FASTQ file... Trying with BAM..."
		#alignment = AlignIO.read(open(FASTA_FILE), "")
		#count = SeqIO.convert(FASTA_FILE, "bam", FASTA_FILE+".fasta" , "fasta")
		#print "Converted %i records" % count
		#print alignment
		#out_file = "%s.fa" % os.path.splitext(FASTA_FILE)[0]
		#with open(out_file, "w") as out_handle:
			# Write records from the BAM file one at a time to the output file.
			# Works lazily as BAM sequences are read so will handle large files.
			#SeqIO.write(bam_to_rec(FASTA_FILE), out_handle, "fasta")	
		for i in bam_to_rec(FASTA_FILE):
			#print "\n\n\n---"
			#print i.seq
			array_seq.append(str(i.seq))
			
		#alignment = AlignIO.read(open(out_file), "fasta")
		#summary_align = AlignInfo.SummaryInfo(alignment)
		#consensus = summary_align.gap_consensus(threshold = 1.0, ambiguous = 'N', consensus_alpha = alphabet, require_multiple = 2)
		#print consensus	
		#print len(consensus)
		#exit(1)	
		pass	
		
	if len(array_seq)==0:
		print "ERROR: empty input"
		exit(1)

	FASTA_FILE="--interactive--"
	print "Creating k-mer in memory (K="+str(K_MER_SIZE)+")..."
	input_fasta_k_mer=convert_to_k_mer_hash(array_seq,min_kmer_repeats)
	print str(len(input_fasta_k_mer))+" K-mers as input"
	#print input_fasta_k_mer
	print "Comparing input with DB..."
	total_hits=0
	
	for input_k_mer in input_fasta_k_mer.keys():
		try:
			#print input_k_mer
			db_line=k_mer_hash[input_k_mer]
			#aux_array_ref=read_var_from_file_line(DB_FILE,db_line)
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
	#max_haplogroup_score_table=aux_array_ref=ast.literal_eval(linecache.getline(DB_FILE,2))
        max_haplogroup_score_table=ast.literal_eval(TEST_content[1])

	#print max_haplogroup_score_table
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
			#print ranking_table[i]
			if result[0]=='':
				result=ranking_table[i]
			else:
				result[0]+=", "+ranking_table[i][0]
				
			score=ranking_table[i][3]
			i+=1
	except IndexError:
		print "ERROR: no result, check input"
		exit(1)
	print result

else:
	print "ERROR: Usage: "+str(sys.argv[0])+" RESULT_DB INPUT_FASTA.fasta"

