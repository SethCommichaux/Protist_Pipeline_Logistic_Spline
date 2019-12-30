import argparse
import os
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-b2b", help="busco2busco mapping file")
parser.add_argument("-e", help="euk_busco.pep file")
parser.add_argument("-b", help="busco_marker_gene_ids.txt file")
parser.add_argument("-f", help="euk functions file")
parser.add_argument("-no", help="ncbi nodes file")
parser.add_argument("-na", help="ncbi names file")
args = parser.parse_args()

lineages = {i.strip().split('\t')[0]:i.strip().split('\t')[3] for i in open(args.f)}
nodes = {i.strip().split('\t')[0]:i.strip().split('\t')[1] for i in open(args.no)}
names = {i.strip().split('\t')[1].upper():i.strip().split('\t')[0] for i in open(args.na)}

print "Built dictionaries successfully!!!"

b2b = {}

for i in open(args.b2b):
	tmp = i.strip().split('\t')
	for h,j in enumerate(tmp):
		if h == 0: b2b[tmp[h]] = tmp[h]
		else:
			b2b[tmp[h]] = tmp[0]

uniref2busco = {i.strip().split('\t')[1]:i.strip().split('\t')[0].split(':')[1] for i in open(args.b)}

print "Built mapping indexes!!!"

for i in SeqIO.parse(args.e,'fasta'):
	id = str(i.id)
	busco = b2b[uniref2busco[id]]
	d = str(i.description)
	s = str(i.seq)
	lineage = map(lambda x: x.strip(),lineages[id].split(';'))
	fgs = {'family':'','genus':'','species':''}
	for j in lineage:
		if nodes[names[j]] in fgs:
			with open(busco+'_'+j.replace(' ','_').replace('/','_')+'.fasta','a') as out:
				out.write(">"+d+"\n"+s+"\n")


for i in os.listdir('.'):
	c = 0
	for j in open(i):
		if j[0] == ">":
			c += 1
	if c == 1:
		os.system('mv '+i+' singletons/')
	elif c > 1:
		os.system('mv '+i+' MSA/')


print "Finished!!!"






## Need to remove marker genes without a least family level taxonomic rank from uniref100 