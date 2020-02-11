from Bio import SeqIO
import argparse
import os

# parse the user provided arguments; all must be present
parser = argparse.ArgumentParser()
parser.add_argument("-p", help="path to directory of msa's")
parser.add_argument("-o", help="output file")
parser.add_argument("-d", help="path to data directory")
args = parser.parse_args()

names = {i.strip().split('\t')[1].upper():i.strip().split('\t')[0] for i in open(args.d+'/names')}

fullLineage = {}

for i in open(args.d+'/fullnamelineage.dmp'):
	tmp = i.upper().strip().split('|')
	id  = tmp[0].strip()
	name = tmp[1].strip()
	lineage = tmp[2].strip() + ' ' + name
	fullLineage[id] = lineage

# collect the length of gene sets per taxa group
gene_lengths = {}

for i in os.listdir(args.p):
	path = args.p+'/'+i
	if path.endswith('.fasta'):
		marker_gene_group_label = i.replace('.fasta','')
		taxaID = marker_gene_group_label.split('_')[1]
		lineage = fullLineage.get(taxaID,0).split(';')
		if lineage == 0: continue
		msa = path+'.msa'
		if os.path.exists(msa) == False:
			for j in SeqIO.parse(path,'fasta'):
				seq_len = len(str(j.seq))
				break
		else:
			for j in SeqIO.parse(msa,'fasta'):
				seq_len = len(str(j.seq))
				break
		for line in lineage:
			line = line.strip()
			line_taxaID = names[line]			
			if line_taxaID not in gene_lengths: gene_lengths[line_taxaID] = seq_len
			else: gene_lengths[line_taxaID] += seq_len
	
with open(args.o,'w') as out:
	for k,v in gene_lengths.items():
		out.write(k+'\t'+str(v)+'\n')
