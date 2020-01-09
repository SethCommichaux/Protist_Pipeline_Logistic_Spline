import argparse
from Bio import SeqIO
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-b", help="busco ids file")
parser.add_argument("-e", help="eukaryota.pep file")
parser.add_argument("-o", help="path name for output fasta file of extracted sequences")
args = parser.parse_args()

id2seq_len = {str(i.id):len(i.seq) for i in SeqIO.parse(args.e,'fasta')}

# extract uniref100 proteins that matched BUSCO marker gene HMM profile; uniref100 proteins are grouped by busco gene
# also record each uniref100 protein's length

busco_candidates = {}

for i in open(args.b):
	tmp = i.strip().split('\t')
	if len(tmp) < 2: 
		print tmp
		continue
	busco = tmp[0]
	uniref = tmp[1]
	prot_len = id2seq_len.get(uniref,0)
	if prot_len == 0: continue
	if busco not in busco_candidates:
		busco_candidates[busco] = {uniref:prot_len}
	else:
		busco_candidates[busco][uniref] = prot_len

print "Finished processing kaiju output file!!!"

# if there is over 2-fold variance in gene lengths for marker gene family
# filter out genes that are over 2 standard deviations away from mean

candidates = {}

for k,v in busco_candidates.items():
	lengths = v.values()
	var = max(lengths)/float(min(lengths))
	sd = np.std(lengths)
	m = np.mean(lengths)
	if var <= 2:
		for z,w in v.items():
			candidates[z] = 0
	else:
		high = m+sd
		low = m-sd
		for z,w in v.items():
			if high > w > low:
				candidates[z] = 0

with open(args.o,'w') as out:
	for i in SeqIO.parse(args.e,'fasta'):
		id = str(i.id)
		if id in candidates:
			out.write(">"+str(i.description)+"\n"+str(i.seq)+"\n")

