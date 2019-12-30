# This script splits uniref100 protein fasta file into 3 files based on taxonomy: protist.pep, fungi.pep, and non_eukaryota.pep
# non_eukaryota.pep does have animal and plant proteins, just not fungi and protist proteins.
# Importantly, if filters out proteins that are unclassified by NCBI taxonomy.

import argparse
from Bio import SeqIO
import random

parser = argparse.ArgumentParser()
parser.add_argument("-g", help="Path to UniRef100 protein fasta file")
parser.add_argument("-t", help="path to NCBI fullnamelineage.dmp taxonomy file")
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

taxaID2_Name_Lineage = build_taxa_dict()


with open('protist.pep','w') as out, open('fungi.pep','w') as out2, open('non_eukaryota.pep','w') as out3:
	for i in SeqIO.parse(args.g,'fasta'):
		d = str(i.description)
		s = str(i.seq)
		taxaID = d.split('TaxID=')[1].split(' ')[0]
		if taxaID not in taxaID2_Name_Lineage: continue
		lineage = taxaID2_Name_Lineage[taxaID]
		if 'EUKARYOTA' in lineage:
			if 'FUNGI' in lineage: out2.write(">"+d+"\n"+s+"\n")
			elif 'METAZOA' in lineage: out3.write(">"+d+"\n"+s+"\n")
			elif 'EMBRYOPHYTA' in lineage: out3.write(">"+d+"\n"+s+"\n")
			elif 'ENVIRONMENTAL SAMPLES' in lineage: continue
			elif 'UNCLASSIFIED EUKARYOTES' in lineage: continue
			else: out.write(">"+d+"\n"+s+"\n")
		else:
			out3.write(">"+d+"\n"+s+"\n")


