# This script splits uniref100 protein fasta file into 3 files based on taxonomy: protist.pep, fungi.pep, and non_eukaryota.pep
# non_eukaryota.pep does have animal and plant proteins, just not fungi and protist proteins.
# Importantly, if filters out proteins that are unclassified by NCBI taxonomy.

import argparse
from Bio import SeqIO
import random

parser = argparse.ArgumentParser()
parser.add_argument("-g", help="Path to UniRef format protein fasta file")
parser.add_argument("-t", help="path to NCBI fullnamelineage.dmp taxonomy file")
parser.add_argument("-no", help="path to nodes taxonomy file")
parser.add_argument("-na", help="path to names taxonomy file")
args = parser.parse_args()


# parse NCBI taxonomy fullnamelineage.dmp file for taxaIDs, taxaNames and taxonomic lineages
def build_taxa_dict():
	taxaID2_Name_Lineage = {}
	for i in open(args.t):
		tmp = i.strip().split('\t|\t')
		taxaID = tmp[0]
		taxaName = tmp[1].upper().strip('\t|\t')
		lineage = tmp[2].upper().strip('\t|\t')
		taxaID2_Name_Lineage[taxaID] = lineage+taxaName
	print "Taxa dictionary built!!!"
	return taxaID2_Name_Lineage

nodes = {i.strip().split('\t')[0]:i.strip().split('\t')[1] for i in open(args.no) if i.strip().split('\t')[1] in ['family','genus','species']}
names = {i.strip().split('\t')[1].upper():i.strip().split('\t')[0] for i in open(args.na) if i.strip().split('\t')[0] in nodes}
taxaID2_Name_Lineage = build_taxa_dict()


with open('eukaryota.pep','w') as out1, open('non_eukaryota.pep','w') as out2:
	for i in SeqIO.parse(args.g,'fasta'):
		d = str(i.description)
		s = str(i.seq)
		taxaID = d.split('TaxID=')[1].split(' ')[0]
		if taxaID not in taxaID2_Name_Lineage: continue
		lineage = taxaID2_Name_Lineage[taxaID]
		use = False
		for group in lineage.split(';'):
			if group.strip() in names:
				use = True
				break
		if use == False: continue
		if 'EUKARYOTA' in lineage:
			if 'FUNGI' in lineage: out1.write(">"+d+"\n"+s+"\n")
			elif 'METAZOA' in lineage: out2.write(">"+d+"\n"+s+"\n")
			elif 'EMBRYOPHYTA' in lineage: out2.write(">"+d+"\n"+s+"\n")
			elif 'ENVIRONMENTAL SAMPLES' in lineage: continue
			elif 'UNCLASSIFIED EUKARYOTES' in lineage: continue
			else: out1.write(">"+d+"\n"+s+"\n")
		else:
			out2.write(">"+d+"\n"+s+"\n")


