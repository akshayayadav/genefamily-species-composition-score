#!/usr/bin/python

import re
import sys
import os
import pandas as pd
import numpy as np
import itertools
from sklearn.metrics.pairwise import cosine_similarity

import argparse

parser = argparse.ArgumentParser(description="Tool for scoring gene families using expected species counts per family")
parser._optionals.title="Arguments"
parser.add_argument('--exp_sp_ct',help="File containing expected gene counts per species", required=True, dest="exp_sp_count_file")
parser.add_argument('--fasta_dir',help="Directory containing family fasta files", required=True, dest="fasta_dir")
parser.add_argument('--out', help="Output file name", required=True, dest="out_file")
args = parser.parse_args()



def read_perfamily_species_count_file(perfamily_species_count_fileName):
	perfamily_species_count_file=open(perfamily_species_count_fileName, "r")
	speciesID_seqcountsArr_dict={}
	for line in perfamily_species_count_file:
		line=line.rstrip()
		speciesID, seqcounts=re.split(r'\:', line)
		speciesID_seqcountsArr_dict[speciesID]=map(int,(re.split(r'\,',seqcounts)))

	
	perfamily_species_count_file.close()
	return(speciesID_seqcountsArr_dict)

def get_species_count_matrix(speciesID_seqcountsArr_dict):
	speciesID_seqcountsExpandedArr_dict={}
	for speciesID in  speciesID_seqcountsArr_dict:
		if(len(speciesID_seqcountsArr_dict[speciesID])==2):
			speciesID_seqcountsExpandedArr_dict[speciesID]=list(range(speciesID_seqcountsArr_dict[speciesID][0], speciesID_seqcountsArr_dict[speciesID][1]+1, 1))
		else:
			speciesID_seqcountsExpandedArr_dict[speciesID]=list(speciesID_seqcountsArr_dict[speciesID])
		
	return(speciesID_seqcountsExpandedArr_dict)

def get_seqcounts_centroid(speciesID_seqcountsExpandedArr_dict):
	speciesID_order_arr=list()
	seqcountsExpanded_arr=[]
	for speciesID in speciesID_seqcountsExpandedArr_dict:
		speciesID_order_arr.append(speciesID)
		seqcountsExpanded_arr.append(speciesID_seqcountsExpandedArr_dict[speciesID])
	
	seqcounts_centroid = calculate_seqcounts_centroid(seqcountsExpanded_arr)
	seqcounts_centroid = list(seqcounts_centroid)
	return([seqcounts_centroid, speciesID_order_arr])
def calculate_seqcounts_centroid(seqcountsExpanded_arr):
	combinations_arr=list()
	for combination in itertools.product(*seqcountsExpanded_arr):
		combinations_arr.append(combination)
	seqcounts_combination_matrix = pd.DataFrame(combinations_arr)
	#seqcounts_combination_matrix=seqcounts_combination_matrix.mean().transpose()
	return(seqcounts_combination_matrix.mean())


def process_fasta_files(family_fasta_dirName, seqcounts_centroid, speciesID_order_arr, output_fileName):
	output_file=open(output_fileName,"w")
	for fasta_fileName in os.listdir(family_fasta_dirName):
		species_seqcount_dict, famsize=get_species_counts_from_family_fasta(family_fasta_dirName+"/"+fasta_fileName)
		species_composition_cosine_score=calculate_species_composition_cosine_score(species_seqcount_dict, seqcounts_centroid, speciesID_order_arr)
		output_file.write(os.path.splitext(os.path.basename(fasta_fileName))[0]+" "+str(famsize)+" " +str(species_composition_cosine_score)+"\n")
		#break

def get_species_counts_from_family_fasta(family_fasta_fileName):
	family_fasta_file=open(family_fasta_fileName,"r")
	species_seqcount_dict={}
	famsize=0
	for line in family_fasta_file:
		line=line.rstrip()
		if not (re.match(r'^\>',line)):
			continue
		speciesID=re.split(r'\.',line)[0]
		speciesID = speciesID[1:]
		if(species_seqcount_dict.has_key(speciesID)):
			species_seqcount_dict[speciesID]+=1
			famsize+=1
		else:
			species_seqcount_dict[speciesID]=1
			famsize+=1
	return([species_seqcount_dict, famsize])

def calculate_species_composition_cosine_score(species_seqcount_dict, seqcounts_centroid, speciesID_order_arr):
	family_seqcount_arr=list()
	for speciesID in speciesID_order_arr:
		if(species_seqcount_dict.has_key(speciesID)):
			family_seqcount_arr.append(species_seqcount_dict[speciesID])
		else:
			family_seqcount_arr.append(0)
	
	family_seqcount_arr = np.array(family_seqcount_arr).reshape(1, -1)
	seqcounts_centroid = np.array(seqcounts_centroid).reshape(1, -1)
	species_composition_cosine_score=cosine_similarity(family_seqcount_arr, seqcounts_centroid)
	species_composition_cosine_score=species_composition_cosine_score[0][0]
	return(species_composition_cosine_score)
#################################################################################################

perfamily_species_count_fileName=args.exp_sp_count_file
family_fasta_dirName=args.fasta_dir
output_fileName=args.out_file

speciesID_seqcountsArr_dict = read_perfamily_species_count_file(perfamily_species_count_fileName)
speciesID_seqcountsExpandedArr_dict = get_species_count_matrix(speciesID_seqcountsArr_dict)
seqcounts_centroid, speciesID_order_arr = get_seqcounts_centroid(speciesID_seqcountsExpandedArr_dict)

process_fasta_files(family_fasta_dirName, seqcounts_centroid, speciesID_order_arr, output_fileName)




