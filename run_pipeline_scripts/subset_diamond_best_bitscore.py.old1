import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-d", type=str,
			help="Input diamond results file (--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen)")
parser.add_argument("-o", type=str,
			help="Desired output path for subsetted Diamond output file")
args = parser.parse_args()


def merge_ranges(ranges):
	ranges = iter(sorted(ranges))
	current_start, current_stop = next(ranges)
	for start, stop in ranges:
		if start > current_stop:
		# Gap between segments: output current segment and start a new one.
			yield current_start, current_stop
			current_start, current_stop = start, stop
		else:
			# Segments adjacent or overlapping: merge.
			current_stop = max(current_stop, stop)
	yield current_start, current_stop


protein_coverage = {}
read2bitscore = {}
results = {}

for i in open(args.d):
	tmp = i.strip().split('\t')
	read,protein,start,end,bitscore = tmp[0],tmp[1],int(tmp[8]),int(tmp[9]),float(tmp[11])
	if read not in read2bitscore:
		read2bitscore[read] = bitscore
		if protein in protein_coverage:
			protein_coverage[protein].append((start,end))
		else:
			protein_coverage[protein] = [(start,end)]

	else:
		if bitscore == read2bitscore[read]:
			if protein in protein_coverage:
				protein_coverage[protein].append((start,end))
			else:
				protein_coverage[protein] = [(start,end)]


for k,v in protein_coverage.items():
	num_reads = len(v)
	x = list(merge_ranges(v))
	coverage = 0
	for z,w in x:
		coverage += w-z
	if float(coverage) >= 0:
		results[k] = True
	elif coverage >= num_reads*30:
			results[k] = True
	else:
		results[k] = False


read2bitscore = {}

with open(args.o,'w') as out:
	for i in open(args.d):
		tmp = i.strip().split('\t')
		read,protein,bitscore = tmp[0],tmp[1],float(tmp[11])
		if read not in read2bitscore:
			read2bitscore[read] = bitscore
			if results.get(protein,False) == True:
				out.write(i)
		else:
			if bitscore == read2bitscore[read]:
				if results.get(protein,False) == True:
					out.write(i)
