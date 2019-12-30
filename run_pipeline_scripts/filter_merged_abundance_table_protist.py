import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-nodes',help='ncbi nodes.dmp')
parser.add_argument('-names',help='ncbi names.dmp')
parser.add_argument('-t',help='input merged abundance table to filter')
args = parser.parse_args()


nodes = {i.strip().split('\t')[0]:i.strip().split('\t')[1] for i in open(args.nodes)}
names = {i.strip().split('\t')[1].upper():i.strip().split('\t')[0] for i in open(args.names)}
merged_table = args.t

rank_results = {'family':[],'genus':[],'species':[]}

for rank in ['family','genus','species']:
	results = {}
	with open(merged_table+'.'+rank,'w') as out:
		for h,i in enumerate(open(merged_table)):
			if h == 0:
				out.write(i)
				samples = i.strip().split('\t')[1:]
			else:
				tmp = i.strip().split('\t')
				lineage = tmp[0].split(';')
				counts = np.array(map(float,tmp[1:]))
				for j in lineage:
					phylogroup = j.strip()
					if phylogroup in names:
						if nodes[names[phylogroup]] == rank.lower():
							rank_results[rank].append(phylogroup)
							if phylogroup not in results:
								results[phylogroup] = counts
							else:
								results[phylogroup] = np.add(results[phylogroup],counts)
	
	
		for k,v in sorted(results.items()):
			out.write(k+'\t'+'\t'.join(map(str,v))+'\n')
		

with open('Ranks_counts.txt','w') as out2:
	for z,w in sorted(rank_results.items()):		
		out2.write(z+'\t'+str(len(set(w)))+'\n')


