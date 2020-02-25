import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str,
			help="Input diamond results file (--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen)")
args = parser.parse_args()

results = {}

with open('multimapped.txt','w') as out, open('singletons.txt','w') as out2:
	for i in open(args.d):
		tmp = i.strip().split('\t')
		read,protein,percent_ident,aln_len,qstart,bitscore = tmp[0],tmp[1],float(tmp[2]),int(tmp[3]),int(tmp[8]),float(tmp[11])
		if aln_len >= 30:
			if read not in results:
				results[read] = {bitscore:[i]}
			elif bitscore in results[read]:
				results[read][bitscore].append(i)

	print results.items()[:10]

	for read,v in results.items():
		for bitscore,proteins in v.items():
			if len(proteins) == 1: out2.write(proteins[0])
			else:
				for x in proteins:
					out.write(x)
