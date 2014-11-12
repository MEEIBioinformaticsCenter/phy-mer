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
from operator import eq


##### DEFAULT VARS
K_MER_SIZE=12
k_mer_hash={}
haplogroup_input={}
haplogroup_fasta={}
number_of_haplogroups=0
DEFAULT_CHANGE_TABLE={'A':'G','G':'A','C':'T','T':'C'}



#### FUNCTIONS

# Array of DNA characters to fasta (string) format
def reconstruct_array_seq(array_seq):
	tmp_string=''
	for i in array_seq:
		tmp_string+=i
	return(list(tmp_string.replace('N',''))) 

# Checking that each haplogroup has at least one k-mer different than all other haplogroups
def check_k_mer_hash(k_mer_hash):
	aux_haplogroup_k_mer={}
	for k_mer in k_mer_hash.keys():
		for haplogroup in k_mer_hash[k_mer].keys():
			if k_mer_hash[k_mer][haplogroup]:
				try:
					aux_haplogroup_k_mer[haplogroup].append(k_mer)
				except KeyError:
					aux_haplogroup_k_mer[haplogroup]=[]
					aux_haplogroup_k_mer[haplogroup].append(k_mer)
					pass
	count=0
	aux_haplogroup_k_mer2=aux_haplogroup_k_mer.copy()
	for haplogroup in aux_haplogroup_k_mer:
		count+=1
		same_k_mer=False
		for haplogroup2 in aux_haplogroup_k_mer2:
			if haplogroup2!=haplogroup:
				if len(aux_haplogroup_k_mer[haplogroup])==len(aux_haplogroup_k_mer[haplogroup2]):
					if all(map(eq,aux_haplogroup_k_mer[haplogroup],aux_haplogroup_k_mer[haplogroup2])):
						same_k_mer=True
						print haplogroup+" is the same as "+haplogroup2
						break

		if same_k_mer:
			break
		del aux_haplogroup_k_mer2[haplogroup]
	return not same_k_mer

# Removing k-mers pressent in all haplogroups
def compress_k_mer_hash(k_mer_hash):
	global number_of_haplogroups
	return_k_mer_hash={}
	for i in k_mer_hash.keys():
		to_return=False
		if len(k_mer_hash[i].keys())<number_of_haplogroups:
			return_k_mer_hash[i]=k_mer_hash[i].copy()
	return return_k_mer_hash

# Function NOT IN USE used for test and debug the code	
def print_fasta_files(haplogroup_fasta):
	for haplogroup_name in haplogroup_fasta.keys():
		wfile=open(haplogroup_name+".FASTA",'w')
        	wfile.write(">gi|"+haplogroup_name+"\n")
		array_seq=reconstruct_array_seq(haplogroup_fasta[haplogroup_name])
		counter=1
		for i in array_seq:
			wfile.write(i)
	        	if counter==70:
				wfile.write("\n")
				counter=0
			counter+=1
		wfile.close()

# Converting a fasta sequence (string) to a dictionary (hash table) of k-mers
def convert_to_k_mer_hash(haplogroup_fasta):
	global K_MER_SIZE
	k_mer_hash={}

	for haplogroup_name in haplogroup_fasta.keys():
		array_seq=reconstruct_array_seq(haplogroup_fasta[haplogroup_name])
		count_base=0
		while count_base<= len(array_seq)-K_MER_SIZE:
			coun_kmer=0
			kmer=''
			while coun_kmer<K_MER_SIZE:
				kmer+=array_seq[count_base+coun_kmer]
				coun_kmer+=1
			try:
				k_mer_hash[kmer][haplogroup_name]=True
			except KeyError:
				k_mer_hash[kmer]={}
				k_mer_hash[kmer][haplogroup_name]=True
				pass
			count_base+=1
	return k_mer_hash


#### MAIN
def main():
	global K_MER_SIZE
	global k_mer_hash
	global haplogroup_input
	global haplogroup_fasta
	global number_of_haplogroups

	if len(sys.argv)<4 or len(sys.argv)>5:
		print "ERROR: Usage: "+str(sys.argv[0])+" REFERENCE_FASTA_FILE.fasta SNPS_HAPLOGROUPS.csv RESULT_DB [min_k-mer_size]"
		exit(1)

	try:
		K_MER_SIZE=sys.argv[4]
	except IndexError:
		pass
		
	FASTA_FILE=sys.argv[1]
	HAPLOGROUP_FILE=sys.argv[2]
	OUTPUT_FILE=sys.argv[3]
	
	print "Reading Haplogroup information..."
	for rline in open(HAPLOGROUP_FILE,'r').readlines():
		rline_splited=rline.replace('\n','').split(',')
		if rline_splited[0]!='Haplogroup':
			haplogroup_input[rline_splited[0]]=[]
			for variant in rline_splited:
				if (variant!=rline_splited[0] and variant!=""):
					haplogroup_input[rline_splited[0]].append(variant)

	haplogroup_array=haplogroup_input.keys()
	number_of_haplogroups=len(haplogroup_array)
	print "Openning fasta file and loading it in memory..."
	

	handle = open(FASTA_FILE, "rU")
	for record in SeqIO.parse(handle, "fasta") :
		array_seq=list(record.seq)
		for aux_haplogroup in haplogroup_array:
			haplogroup_fasta[aux_haplogroup]=list(array_seq)
		break
	handle.close()

	print "Adding mutations to each Haplogroup..."
	for aux_haplogroup in haplogroup_array:
		for aux_mutation in haplogroup_input[aux_haplogroup]:
			try: ### Is a number, so standar change
				pos=int(aux_mutation)-1
				try:
					haplogroup_fasta[aux_haplogroup][pos]=DEFAULT_CHANGE_TABLE[haplogroup_fasta[aux_haplogroup][pos]]
				except KeyError:
					print "Warning!!! impossible to apply:"
					print "aux_haplogroup: "+str(aux_haplogroup)
					print "pos: "+str(pos)
					print "Ref: "+str(haplogroup_fasta[aux_haplogroup][pos])
					#exit(1)
					pass
				#print haplogroup_fasta[aux_haplogroup][pos]
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
						haplogroup_fasta[aux_haplogroup][pos]+=insertion
				elif aux_mutation.find('d')!=-1: #### there is a 'd' so it's a deletion
					if aux_mutation.find('-')!=-1: ### Range	
						str_start,str_end=aux_mutation.replace('d','').split('-')
						start=int(str_start)-1
						end=int(str_end)-1
						for pos in range(start,end+1):
							haplogroup_fasta[aux_haplogroup][pos]=''
					else : ### only one base
						try:
							pos=int(aux_mutation.replace('d',''))-1
							haplogroup_fasta[aux_haplogroup][pos]=''
						except ValueError: ## Deleting an specific letter
							ref_to_delete=aux_mutation[0]
							pos=int(aux_mutation.replace('d','').replace(ref_to_delete,''))-1
							if haplogroup_fasta[aux_haplogroup][pos]!=ref_to_delete:
								print "WARNING!!!!!! Removing "+haplogroup_fasta[aux_haplogroup][pos]+" instead of "+ref_to_delete+" at "+aux_mutation
                                                        haplogroup_fasta[aux_haplogroup][pos]=''
							pass
				else: ###### Number with letter!
					snp=aux_mutation[-1]
					pos=int(aux_mutation[0:-1])-1
					haplogroup_fasta[aux_haplogroup][pos]=snp
				pass
	k_mer_hash={}
	
	valid_k_mer=False
	while not valid_k_mer:
		k_mer_hash.clear()
		print "Creating k-mer DB in memory with k="+str(K_MER_SIZE)
		k_mer_hash=convert_to_k_mer_hash(haplogroup_fasta)
		print len(k_mer_hash.keys())
		print "Compressing DB in memory..."
		k_mer_hash=compress_k_mer_hash(k_mer_hash)
		print len(k_mer_hash.keys())
		print "Checking that all Haplogroups are unique..."
		valid_k_mer=check_k_mer_hash(k_mer_hash)
		if valid_k_mer:
			print "Preparing data to be stored in a plain text file..."
			output_kmer_hash={}
			output_ref_hash={}
			output_ref_index=[]
			output_ref_index_opposite_hash={}
			lines_array=[]

			count=3 ##### first 3 lines with index
			for kmer in k_mer_hash.keys():
				count+=1
				output_kmer_hash[kmer]=count
				lines_array.append(k_mer_hash[kmer].keys())
				for haplogroup in k_mer_hash[kmer].keys():
					try:
						output_ref_hash[haplogroup]+=1
					except KeyError:
						output_ref_hash[haplogroup]=1
			#### Changing haplogroups with numbers!!!!!
			i=0
			while i<len(lines_array):
				new_line=[]
				for element in lines_array[i]:
					try:
						new_line.append(output_ref_index_opposite_hash[element])
					except KeyError:
						output_ref_index.append(element)
						output_ref_index_opposite_hash[element]=len(output_ref_index)-1
						new_line.append(output_ref_index_opposite_hash[element])
						
						pass
				lines_array[i]=list(new_line)
				i+=1
			print "Saving data in: "+OUTPUT_FILE
	
			f = open(OUTPUT_FILE, 'w')
			f.write(str(output_kmer_hash)+"\n")
			f.write(str(output_ref_hash)+"\n")
			f.write(str(output_ref_index)+"\n")
			for line in lines_array:
				f.write(str(line).replace('[','').replace(']','')+"\n")
			f.close	
		else:
			K_MER_SIZE+=1
			print "K-mer to small, changing to "+str(K_MER_SIZE)
			


##############
if __name__ == "__main__":
    main()


