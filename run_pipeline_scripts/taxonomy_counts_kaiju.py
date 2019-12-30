import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f", help="euk_functions.txt file")
parser.add_argument("-k", help="kaiju output file")
parser.add_argument("-o", help="path name for output file")
args = parser.parse_args()

uniref2taxaLineage = {i.strip().split('\t')[0]:i.strip().split('\t')[3] for i in open(args.f)}
results = {}

for i in open(args.k):
	tmp = i.strip().split('\t')
	unirefs = tmp[3].split(',')
	for j in unirefs:
		if j in uniref2taxaLineage:
			lineage = uniref2taxaLineage[j]
			if lineage not in results:
				results[lineage] = [1,[j]]
			else:
				results[lineage][0] += 1
				results[lineage][1].append(j)


for k,v in results.items():
	results[k] = [v[0],len(set(v[1]))]


with open(args.o,'w') as out:
	for k,v in sorted(results.items(),key=lambda x:x[1][1]):
		read_count = v[0]
		protein_count = v[1]
		if protein_count >= 2:
			out.write(k+'\t'+str(v[0])+'\t'+str(v[1])+'\n')

