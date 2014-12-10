#!/usr/bin/env python

#    Phy-Mer.py
#    Copyright (C) 2014  Daniel Navarro-Gomez (Daniel_navarro@meei.harvard.edu)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.







import sys
from operator import itemgetter
import ast
import os
import pysam
import getopt
from Bio import SeqIO, Seq, SeqRecord

#### GLOBAL VARS:
K_MER_SIZE=0
REF_INDEX_ARRAY=''
PRINT_RANKING=False
min_kmer_repeats_bam=10


IUPAC_ambiguity_dic={}
IUPAC_ambiguity_dic['U']=['T']
IUPAC_ambiguity_dic['M']=['A','C']
IUPAC_ambiguity_dic['R']=['A','G']
IUPAC_ambiguity_dic['W']=['A','T']
IUPAC_ambiguity_dic['S']=['C','G']
IUPAC_ambiguity_dic['Y']=['C','T']
IUPAC_ambiguity_dic['K']=['G','T']
IUPAC_ambiguity_dic['V']=['A','C','G']
IUPAC_ambiguity_dic['H']=['A','C','T']
IUPAC_ambiguity_dic['D']=['A','G','T']
IUPAC_ambiguity_dic['B']=['C','G','T']
#IUPAC_ambiguity_dic['N']=['A','C','G','T']  ##You should include this one only in case that you do NOT use N to fill up non-sequenced regions.

#### FUNCTIONS
# convert a string in a array of strings using all posibilities for IUPAC
def decompress_IUPAC_kmer(kmer):
	global IUPAC_ambiguity_dic
	result_array=[kmer]
	change=True
	while change:
		change=False
		aux_kmer_counter=0
		while aux_kmer_counter < len(result_array):
			for iupac_code in IUPAC_ambiguity_dic.keys():
				if result_array[aux_kmer_counter].find(iupac_code)!=-1: 
					change=True
					first_iupac=True
					for iupac_code_content in IUPAC_ambiguity_dic[iupac_code]:
						if first_iupac:
							first_iupac=False
							aux_code=result_array[aux_kmer_counter]
							result_array[aux_kmer_counter]=aux_code.replace(iupac_code,iupac_code_content,1)
						else:
							result_array.append(aux_code.replace(iupac_code,iupac_code_content,1))

			aux_kmer_counter+=1
	return(result_array)
					

# Converting a fasta sequence (string) to a dictionary (hash table) of k-mers
def convert_to_k_mer_hash(array_seq,min_repeats=1):
	global IUPAC_ambiguity_dic
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
			
			## IUPAC modification
			IUPAC_used=False
			for iupac_code in IUPAC_ambiguity_dic.keys():
				if kmer.find(iupac_code)!=-1:
					IUPAC_used=True
			if IUPAC_used:
				for aux_kmer in decompress_IUPAC_kmer(kmer):
					try:	
						return_k_mer[aux_kmer]+=1
					except KeyError:
						return_k_mer[aux_kmer]=1
			##

			else:
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

# Reads a snp_def_file in memory to annotate results
def read_snp_def_file(csv_file):
	haplogroup_input={}
        for rline in open(csv_file,'r').readlines():
                rline_splited=rline.replace('\n','').split(',')
                if rline_splited[0]!='Haplogroup':
                        haplogroup_input[rline_splited[0]]=[]
                        for variant in rline_splited:
                                if (variant!=rline_splited[0] and variant!=""):
                                        haplogroup_input[rline_splited[0]].append(variant)
	return(haplogroup_input)


# Help function
def print_help():
	print "Usage: "+str(sys.argv[0])+" [--verbose] [--print-ranking] [--def-snps=haplogroup_def_motifs.csv] [--min-DoC="+str(min_kmer_repeats_bam)+"] Library.txt INPUT_1 [INPUT_2 ... INPUT_X]"
	print "Novel mitochondrial genome haplogroup defining algorithm using a k-mer approach."
	print ""
	print "Optinal arguments."
	print "  --help                     Show this help."
	print "  --verbose                  Print step by step process."
	print "  --print-ranking            Print top 5 results instead of the best match."
	print "  --def-snp=file.csv         Add Haplogroup defining snps to top matches"
	print "                             based in file.csv (Build_16_-_rCRS-based_haplogroup_motifs.csv"
	print "                             in resources folder)."
	print "  --min-DoC=10               Only apply to BAM inputs. Minimal number of occurences of a K-mer"
	print " 			    to be consider."
	print ""


#### MAIN
def main():
	verbose=False
	DEF_SNP=''
	global PRINT_RANKING
	global min_kmer_repeats_bam
	try:
		opts, args = getopt.getopt(sys.argv[1:], '', ['verbose','print-ranking','help','def-snps=','min-DoC='])
	except getopt.GetoptError:
		print "ERROR: Usage: "+str(sys.argv[0])+" [--verbose] [--print-ranking] [--def-snps=haplogroup_def_motifs.csv] [--min-DoC="+str(min_kmer_repeats_bam)+"] Library.txt INPUT_1 [INPUT_2 ... INPUT_X]"
		print "Use --help for more informtion."
                exit(1)
	
	for o,p in opts:
		if o in ['--min-DoC']:
			min_kmer_repeats_bam=int(p)
		if o in ['--verbose']:
			verbose=True
		if o in ['--print-ranking']:
			PRINT_RANKING=True
		if o in ['--def-snps']:
			DEF_SNP=p
		if o in ['--help']:
			print_help()
			exit(0)

	global REF_INDEX_ARRAY
	global K_MER_SIZE
	if len(args)<2:
		print "ERROR: Usage: "+str(sys.argv[0])+" [--verbose] [--print-ranking] [--def-snps=haplogroup_def_motifs.csv] [--min-DoC="+str(min_kmer_repeats_bam)+"] DataBase.txt INPUT_1 [INPUT_2 ... INPUT_X]"
		print "Use --help for more informtion."
		exit(1)
	if verbose and DEF_SNP!='':
                print "Opening haplogroup defining snps..."
	if DEF_SNP!='':
		haplogroup_snp_dict=read_snp_def_file(DEF_SNP)
	else:
		haplogroup_snp_dict={}
	min_kmer_repeats=1
	FASTA_FILES=args[1:]
	Library_FILE=args[0]
	if verbose:	
		print "Opening Library and Checking k-mer size..."	
		print "READING WHOLE FILE IN MEMORY"
	with open(Library_FILE) as f:
		TEST_content = f.readlines()
	if verbose:
		print "DONE"
	
	
	REF_INDEX_ARRAY=ast.literal_eval(TEST_content[2])
		
	k_mer_hash=ast.literal_eval(TEST_content[0])
	K_MER_SIZE=len(k_mer_hash.keys()[1])
	if verbose:
		print "K-mer="+str(K_MER_SIZE)
	
	for FASTA_FILE in FASTA_FILES:


		if verbose:		
			print "Opening input file and loading it in memory..."
		try:
			handle = open(FASTA_FILE, "rU")
		except IOError:
			print FASTA_FILE+" doesn't exist..."
			exit(1)

		if verbose:	
			print "Processing "+FASTA_FILE	
		total_hits=0
		result_haplogroup_ranking={}
		array_seq=[]
		
		if FASTA_FILE.split('.')[-1].lower()=="fasta" or FASTA_FILE.split('.')[-1].lower()=="fa": 
			if verbose:	
				print "Detected FASTA format"	

			try:
				for record in SeqIO.parse(handle, "fasta") :
					array_seq.append(list(record.seq))
				handle.close()
				if len(array_seq)==0:
					exit(1)
			except:
				print "ERROR: reading FASTA file, please check input"
				exit(1)
				
		elif FASTA_FILE.split('.')[-1].lower()=="fastq": 
			if verbose:	
				print "Detected FASTQ format"	
			try:
				handle = open(FASTA_FILE, "rU")
			except IOError:
				print FASTA_FILE+" doesn't exist..."
				exit(1)	
			array_seq=[]
			for record in SeqIO.parse(handle, "fastq") :
				array_seq.append(list(record.seq))
			handle.close()


		elif FASTA_FILE.split('.')[-1].lower()=="bam": 
			if verbose:	
				print "Detected BAM format"	
			try:
				min_kmer_repeats=min_kmer_repeats_bam
				array_seq=[]
				for i in bam_to_rec(FASTA_FILE):
					array_seq.append(str(i.seq))
			except:
				print "ERROR: reading BAM file, please check input"
				exit(1)
		
		else: 
			print "ERROR: Input file not compatible: fasta, fastq or bam"
			exit(1)

		
		if len(array_seq)==0:
			print "ERROR: empty input"
			exit(1)
		
		if verbose:
			print "Creating k-mer in memory (K="+str(K_MER_SIZE)+")..."
		input_fasta_k_mer=convert_to_k_mer_hash(array_seq,min_kmer_repeats)
		if verbose:
			print str(len(input_fasta_k_mer))+" K-mers as input"
			print "Comparing input with Library..."
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
		if verbose:
			print "Creating score table..."
	        max_haplogroup_score_table=ast.literal_eval(TEST_content[1])
	
		ranking_table=[]
		for haplogroup in result_haplogroup_ranking.keys():
			ranking_table.append([haplogroup,float(result_haplogroup_ranking[haplogroup])/float(max_haplogroup_score_table[haplogroup]),float(result_haplogroup_ranking[haplogroup])/float(total_hits),((float(result_haplogroup_ranking[haplogroup])/float(max_haplogroup_score_table[haplogroup]))+(float(result_haplogroup_ranking[haplogroup])/float(total_hits)))/2])
		ranking_table.sort(key=itemgetter(3),reverse=True)
	
		i=0
		score=0.00
		result=['',0.00,0.00,0.00]
		try:
			if not PRINT_RANKING:
				snps=''
				try:
					while i<10:  ## we are printing only top scores
						if score>ranking_table[i][3] :
							break
						if result[0]=='':
							result=ranking_table[i]
							snps="\t["+str(haplogroup_snp_dict[ranking_table[i][0]])
						else:
							result[0]+=", "+ranking_table[i][0]
							snps+=', '+str(haplogroup_snp_dict[ranking_table[i][0]])
							
						score=ranking_table[i][3]
						i+=1
					snps+=']'
				except KeyError:
					score=0.00
					result=['',0.00,0.00,0.00]
					while i<10:  ## we are printing only top scores
						if score>ranking_table[i][3] :
							break
						if result[0]=='':
							result=ranking_table[i]
						else:
							result[0]+=", "+ranking_table[i][0]
							
						score=ranking_table[i][3]
						i+=1
					pass

				if verbose:
					print str(FASTA_FILE)+"\t"+str(result)+snps
				else:
					print str(FASTA_FILE)+"\t"+str(result[0])+"\t"+str(result[3])+snps
					
			else:
				print str(FASTA_FILE)
				while i<5:
				#while i<len(ranking_table):
					if verbose:
						try:
							print str(ranking_table[i])+"\t"+str(haplogroup_snp_dict[ranking_table[i][0]])
						except KeyError:
							print str(ranking_table[i])
							pass
					else:
						try:
							print str(ranking_table[i][0])+"\t"+str(ranking_table[i][3])+"\t"+str(haplogroup_snp_dict[ranking_table[i][0]])
						except KeyError:
							print str(ranking_table[i][0])+"\t"+str(ranking_table[i][3])
							pass
					i+=1
					
		except IndexError:
			print str(FASTA_FILE)+" ERROR: no result, check input file"
			#exit(1)
	
##############
if __name__ == "__main__":
	main()


