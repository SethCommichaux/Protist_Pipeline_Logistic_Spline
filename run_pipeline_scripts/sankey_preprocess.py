import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-nodes", type=str,
			help="NCBI taxonomy nodes.dmp file")
parser.add_argument("-lineage", type=str,
			help="NCBI taxonomy fullnamelineage.dmp file")
parser.add_argument("-results", type=str,
			help="Taxonomic report file")
parser.add_argument("-o", type=str,
			help="Desired path to output directory")
args = parser.parse_args()


nodes,links,c = {"cellular organisms":0},{},0

name2taxid = {i.split("\t|\t")[1].upper():i.split("\t|\t")[0] for i in open(args.lineage)} #lineage
taxid2taxaGroup = {i.split("\t|\t")[0]:i.split("\t|\t")[2] for i in open(args.nodes)} #nodes
ranks = {"superkingdom":0,
			"kingdom":0,
			"phylum":0,
			"class":0,
			"family":0,
			"genus":0,
			"species":0}

for h,i in enumerate(open(args.results)):
	if h == 0: continue
	zz = i.strip().split('\t')
	lineage = map(lambda x: x.strip(), zz[0].split(';'))
	switch = 0
	for k,j in enumerate(lineage):
		if j == '': continue
		j = j.upper()
		taxid = name2taxid[j]
		rank = taxid2taxaGroup[taxid].lower()
		if taxid == "131567":
			rank = "cellular organisms"
		if rank in ranks:
			switch += 1
			if j not in nodes:
				c += 1
				nodes[j] = c
			if switch > 1:
				link = str(nodes[tmp])+'x'+str(nodes[j])
				if link in links:
					links[link] += 1
				else:
					links[link] = 1
			elif switch == 1:
				link = str(0)+'x'+str(nodes[j])
				if link in links:
					links[link] += 1
				else:
					links[link] = 1

			tmp = j

with open(args.o+'/nodes','w') as out:
	for k,v in sorted(nodes.items(),key=lambda x:x[1]):
		out.write(str(k)+"\t"+str(v)+"\n")

with open(args.o+'/links','w') as out:
	for k,v in links.items():
		out.write(str(k.replace('x','\t'))+"\t"+str(v)+"\n")

