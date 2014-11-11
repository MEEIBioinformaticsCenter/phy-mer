#!/usr/bin/env python

#    build_Phy-Mer_DB.py
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



#### IMPORTS
import sys,getopt,commands,string
from Bio import SeqIO
#from operator import eq


##### DEFAULT VARS
haplogroup_input={}
haplogroup_fasta={}
DEFAULT_CHANGE_TABLE={'A':'G','G':'A','C':'T','T':'C'}



#### FUNCTIONS

# Array of DNA characters to fasta (string) format
def reconstruct_array_seq(array_seq):
	tmp_string=''
	for i in array_seq:
		tmp_string+=i
	return(tmp_string) 

#### MAIN
def main():
	global haplogroup_fasta

	if len(sys.argv)!=4:
		print "ERROR: Usage: "+str(sys.argv[0])+" REFERENCE_FASTA_FILE.fasta SNP.csv RESULT.fasta"
		exit(1)

	FASTA_FILE=sys.argv[1]
	SNP_FILE=sys.argv[2]
	OUTPUT_FILE=sys.argv[3]
	input_snps=[]	

	print "Reading SNP information..."
	for rline in open(SNP_FILE,'r').readlines():
		rline_splited=rline.replace('\n','').split(',')
		haplogroup_input[rline_splited[0]]=[]
		for variant in rline_splited:
			if (variant!=""):
				input_snps.append(variant)

	print "Openning fasta file and loading it in memory..."
	handle = open(FASTA_FILE, "rU")
	for record in SeqIO.parse(handle, "fasta") :
		array_ref_seq=list(record.seq)
		break
	handle.close()
	
	array_result_seq=list(array_ref_seq)

	print "Adding mutations to reference fasta.."
	for aux_mutation in input_snps:
		try: ### Is a number, so standar change
			pos=int(aux_mutation)-1
			try:
				array_result_seq[pos]=DEFAULT_CHANGE_TABLE[array_ref_seq[pos]]
			except KeyError:
				print "Warning!!! impossible to apply:"
				print "aux_haplogroup: "+str(aux_haplogroup)
				print "pos: "+str(pos)
				print "Ref: "+str(array_result_seq[pos])
				#exit(1)
				pass
			#print array_result_seq[pos]
		except ValueError:
			if aux_mutation.find('.')!=-1: #### there is a '.' so it's an insertion
				pos=int(aux_mutation.split('.')[0])-1
				times_string=''
				for i in list(aux_mutation.split('.')[1]):
					try:
						times_string+=str(int(i))
					except ValueError:
						break
				try:
					times=int(times_string)
				except ValueError:
					####### PATCH FOR X!!!!!!!!!!!
					times_string='X'
					times=1
					pass
					####### PATCH FOR X!!!!!!!!!!!
					
				insertion=aux_mutation.split('.')[1].replace(times_string,'')
				for i in range(0,times):
					array_result_seq[pos]+=insertion
			elif aux_mutation.find('d')!=-1: #### there is a 'd' so it's a deletion
				if aux_mutation.find('-')!=-1: ### Range	
					str_start,str_end=aux_mutation.replace('d','').split('-')
					start=int(str_start)-1
					end=int(str_end)-1
					for pos in range(start,end+1):
						array_result_seq[pos]=''
				else : ### only one base
					try:
						pos=int(aux_mutation.replace('d',''))-1
						array_result_seq[pos]=''
					except ValueError: ## Deleting an specific letter
						ref_to_delete=aux_mutation[0]
						pos=int(aux_mutation.replace('d','').replace(ref_to_delete,''))-1
						if array_result_seq[pos]!=ref_to_delete:
							print "WARNING!!!!!! Removing "+array_result_seq[pos]+" instead of "+ref_to_delete+" at "+aux_mutation
						array_result_seq[pos]=''
						pass
			else: ###### Number with letter!
				snp=aux_mutation[-1]
				pos=int(aux_mutation[0:-1])-1
				array_result_seq[pos]=snp
			pass
	
	print "Saving data in: "+OUTPUT_FILE
	f = open(OUTPUT_FILE, 'w')
	f.write(">gi|fasta_from:"+SNP_FILE+"\n")
	chunk_size=70
	result_seq=reconstruct_array_seq(array_result_seq)
	chunk_size
	while i<len(result_seq):
		f.write(result_seq[i:i+chunk_size]+"\n")
		i+=chunk_size
	f.close	

##############
if __name__ == "__main__":
    main()


