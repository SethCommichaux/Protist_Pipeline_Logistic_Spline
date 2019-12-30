import sys
import os

taxa = {}
directory = sys.argv[1]
samples = []

for i in sorted(os.listdir(directory)):
	if i.endswith('.fastq'):
		f = i.replace('.fastq','')
		for j in open(f+'/Taxonomic_Report.txt'):
			tmp = j.strip().split('\t')
			phylogroup = tmp[0]
			taxa[phylogroup] = []



for i in sorted(os.listdir(directory)):
	counts = {}
	if i.endswith('.fastq'):
		f = i.replace('.fastq','')
		samples.append(f)
		for j in open(f+'/Taxonomic_Report.txt'):
			tmp = j.strip().split('\t')
			phylogroup = tmp[0]
			count = float(tmp[1])
			abund = float(tmp[2])
			if count > 1:
				counts[phylogroup] = abund

		for k,v in taxa.items():
			if k in counts: taxa[k].append(counts[k])
			else: taxa[k].append(0)


with open('merged_abundance_table.txt','w') as out:
	out.write('Taxa')
	for j in samples: out.write('\t'+j)
	out.write('\n')
	for z,w in sorted(taxa.items()):
		if w.count(0) == len(w): continue
		out.write(z+'\t'+'\t'.join(map(str,w))+'\n')

