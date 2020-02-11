from Bio import SeqIO
import argparse
import os

# parse the user provided arguments; all must be present
parser = argparse.ArgumentParser()
parser.add_argument("-nodes", help="ncbi nodes.dmp file")
parser.add_argument("-p", help="path to directory of msa's")
parser.add_argument("-o", help="output file")
args = parser.parse_args()

def process_fasta_file(results,msa,taxa_rank,marker_gene_group_label):
	for j in SeqIO.parse(msa,'fasta'):
		id = str(j.id)
		seq = str(j.seq)
		if id not in results:
			results[id] = {'F':[[],[],0],'G':[[],[],0],'S':[[],[],0]}
		results[id][taxa_rank][2] = marker_gene_group_label
		c = 0
		for g,b in enumerate(seq):
			if b == '-':
				if c > 0: c+=1
				else:
					results[id][taxa_rank][0].append(str(g))
					c +=1
			else:
				if c != 0:
					results[id][taxa_rank][1].append(str(c))
					c = 0
		if c != 0: results[id][taxa_rank][1].append(str(c))
	return results

# create dictionaries for ncbi taxonomy nodes files
# nodes = map taxid to taxonomic rank
nodes = {i.strip().split('\t')[0]:i.strip().split('\t')[1] for i in open(args.nodes)}

results = {}

# iterate through directory with marker gene fasta and multiple sequence alignment files
for h,i in enumerate(os.listdir(args.p)):
	path = args.p+'/'+i
	if path.endswith('.fasta'):
		marker_gene_group_label = i.replace('.fasta','')
		taxaID = marker_gene_group_label.split('_')[1]
		rank = nodes.get(taxaID,0)
		taxa_rank = {'family':'F','genus':'G','species':'S'}.get(rank,0)
		if taxa_rank == 0: continue
		msa = path+'.msa'
		if os.path.exists(msa) == False:
			results = process_fasta_file(results,path,taxa_rank,marker_gene_group_label)
		else:
			results = process_fasta_file(results,msa,taxa_rank,marker_gene_group_label)

with open(args.o,'a') as out:
	for uniref,v in results.items():
		out.write(uniref)
		for rank,gaps in sorted(v.items()):
			marker_gene_group_label = gaps[2]
			if gaps[0] == []: gap_starts,gap_lengths = '0','0'
			else:
				gap_starts = ','.join(gaps[0])
				gap_lengths = ','.join(gaps[1])
			out.write('\t'+str(marker_gene_group_label)+'\t'+str(gap_starts)+'\t'+str(gap_lengths))
		out.write('\n')

