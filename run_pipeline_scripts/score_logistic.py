# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen

import argparse
import collections
import numpy as np
from sklearn.linear_model import LogisticRegression

def LCA(taxa):
	for i,j in enumerate(taxa):
		c = 0
		j = j.split(';')
		if i == 0:
			new_lineage = j
		else:
			for i in range(min(len(j),len(new_lineage))):
				if j[i] == new_lineage[i]:
					c += 1
			new_lineage = new_lineage[:c]
	new_lineage = ';'.join(new_lineage)
	return new_lineage


parser = argparse.ArgumentParser()
parser.add_argument("-d", help="diamond output file")
parser.add_argument("-f", help="euk functions file")
parser.add_argument("-s", help="slopes file for marker genes")
parser.add_argument("-no", help="ncbi nodes file")
parser.add_argument("-na", help="ncbi names file")
parser.add_argument("-o", help="path name for output file")
args = parser.parse_args()

lineages = {i.strip().split('\t')[0]:i.strip().split('\t')[3] for i in open(args.f)}
nodes = {i.strip().split('\t')[0]:i.strip().split('\t')[1] for i in open(args.no)}
names = {i.strip().split('\t')[1].upper():i.strip().split('\t')[0] for i in open(args.na)}

first_pass_collect_protein_start = {}

for i in open(args.d):
	tmp = i.strip().split('\t')
	uniref = tmp[1]
	start = int(tmp[8])
	if uniref not in first_pass_collect_protein_start:
		first_pass_collect_protein_start[uniref] = {start}
	else:
		first_pass_collect_protein_start[uniref].update({start})

print "\nExtracted diamond output information!!!"


recruit_protein_data = {}

for i in open(args.s):
	tmp = i.strip().split('\t')
	marker_gene = tmp[0]
	if marker_gene in first_pass_collect_protein_start:
		start = int(tmp[1])
		slopes = map(float,tmp[2].split(','))
		labels = tmp[3].split(',')
		for j in first_pass_collect_protein_start[marker_gene]:
			if abs(start - j) <= 40:
				if marker_gene+'_'+str(j) not in recruit_protein_data:
					recruit_protein_data[marker_gene+'_'+str(j)] = [slopes,labels]
				break

print "\nExtracted marker gene classifier information!!!"

read_taxonomy = {}

with open('read_classified.txt','w') as out, open('candidate_reads.fasta','w') as out2:
	for i in open(args.d):
		tmp = i.strip().split('\t')
		uniref = tmp[1]
		start = int(tmp[8])
		lookup = uniref+'_'+str(start)
		if lookup in recruit_protein_data:
			read = tmp[0]
			aln_read_seq = tmp[12]
			aln_len = float(tmp[3])
			bitscore = float(tmp[11])
			slope = bitscore/aln_len
			x_train = recruit_protein_data[lookup][0]
			y_train = recruit_protein_data[lookup][1]
			if min(len(set(x_train)),len(set(y_train))) < 2:
				continue
			x_train = np.array(x_train).reshape(-1, 1)
			y_train = np.array(y_train)
			lr = LogisticRegression().fit(x_train, y_train)
			coef = float(lr.coef_[0][0])
			y_int = float(lr.intercept_[0])
			probabilities = list(lr.predict_proba(slope)[0])
			lineage = lineages[uniref]
			labels = list(lr.classes_)
			out.write(read+'\t'+uniref+'\t'+','.join(labels)+'\t'+','.join(map(str,probabilities))+'\t'+lineage+'\n')
			tmp = {labels[i]:probabilities[i] for i in range(len(labels))}
			if tmp.get('N',0) >= 0.5: continue
			else:
				try: del tmp['N']
				except: print "No negative class!"
				out2.write(">"+read+"\n"+aln_read_seq+"\n")
				for k,v in sorted(tmp.items(),key=lambda x:-x[1]):
					taxa = ''
					for i in lineage.split(';'): # this extends lineage to correct rank
						i = i.strip()
						if {'species':'S','genus':'G','family':'F'}.get(nodes[names[i]],0)!= k:
							taxa += i+';'
						else:
							taxa += i+';'
							break
					x = 0
					for i in 'FGS':
						if i != k:
							x += tmp.get(i,0)
						else:
							x += tmp[i]
							break
					if read not in read_taxonomy:
						read_taxonomy[read] = [[taxa],[v]]
					else:
						read_taxonomy[read][0].append(taxa)
						read_taxonomy[read][1].append(v)
					break

print "\nFinished classification of reads!!!"
print "\nNumber of positively classified reads: ",len(read_taxonomy)

taxa_counts = {}

with open(args.o,'w') as out:
	for k,v in read_taxonomy.items():
		if v[1].count(max(v[1])) == 1:
			lineage = v[0][v[1].index(max(v[1]))]
			taxa = ''
			for j in lineage.split(';')[:-1]:
				taxa += j+';'
				if taxa not in taxa_counts:
					taxa_counts[taxa] = max(v[1])
				else:
					taxa_counts[taxa] += max(v[1])
		else:
			x = [v[0][i] for i in range(len(v[0])) if v[1][i] == max(v[1])]
			lineage = LCA(x)
			taxa = ''
			for j in lineage.split(';'):
				taxa += j+';'
				if taxa not in taxa_counts:
					taxa_counts[taxa] = max(v[1])
				else:
					taxa_counts[taxa] += max(v[1])

	for k,v in sorted(taxa_counts.items()):
		out.write(k+'\t'+str(v)+'\n')







