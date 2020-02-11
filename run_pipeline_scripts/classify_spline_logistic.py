import os
import math
import argparse
import itertools

# transforms the local start position of the read mapped to the marker gene
# to the global start position in the MSA of the taxonomic group for the marker gene
def adjust_start(start,gaps):
	for gap_start,gap_len in zip(*gaps):
		if start >= gap_start:
			start += gap_len
	return start

# takes the coefficients of the spline polynomial and logistic regression
# models, as well as the classification power of the given region, and returns the
# probability that the aligned read comes from that taxonomic group, as well as the
# classification power of the region
def classify(dist_from_knot,qslope,A,B,C,D,LR_slope,LR_y_int,Classification_power):
	if LR_slope <= 0:
		return 0,Classification_power
	decision_threshold = ((A*dist_from_knot)**3)+((B*dist_from_knot)**2)+(C*dist_from_knot)+D
	if LR_slope == 'NA':
		if qslope >= decision_threshold:
			return 1,Classification_power
		else:
			return 0,Classification_power
	else:
		dist_from_threshold_boundary = qslope - decision_threshold
		probability = 1/(1+(math.exp(-1*(LR_y_int+LR_slope*dist_from_threshold_boundary))))
		return probability,Classification_power

# returns the lowest common ancestor of a list of NCBI taxonomic lineages
def LCA(lineages):
	lineages = map(lambda x: x.split(';'),lineages)
	nr_lineages = []
	for a, b in itertools.combinations(lineages, 2): # if there is best hit tie between different ranks of same lineage, prefer more specific rank
		if len(set(a) | set(b)) == max(a,b,key=len):
			nr_lineages.append(max(a,b,key=len))
		else:
			nr_lineages.append(a)
			nr_lineages.append(b)
	new_lineage = [i[0].strip() for i in zip(*nr_lineages) if len(set(i)) == 1 if '' not in i]
	return new_lineage

# parse user input argumets
parser = argparse.ArgumentParser()
parser.add_argument("-d", help="diamond output file")
parser.add_argument("-path", help="path to data directory i.e. /usr/local/protist_pipeline/data/")
parser.add_argument("-th", help="Probability cutoff for reporting results")
args = parser.parse_args()

# build dictionaries for looking up taxonomic lineages, ID's and ranks
nodes = {i.strip().split('\t')[0]:i.strip().split('\t')[1] for i in open(args.path+'/nodes')}
names = {i.strip().split('\t')[1].upper():i.strip().split('\t')[0] for i in open(args.path+'/names')}
fullLineage = {}

for i in open(args.path+'/fullnamelineage.dmp'):
	tmp = i.strip().split('\t|\t')
	id  = tmp[0].strip()
	name = tmp[1].strip().upper()
	lineage = tmp[2].strip('\t|').upper().strip() + ' ' + name
	fullLineage[id] = lineage

# build dictionaries for looking up taxonomic lineages by uniref ID and protein lengths by marker gene group ID
uniref2lineage = {}
busco2protlen = {}

with open(args.path+'/euk_functions.txt') as f:
	next(f)
	for i in f:
		tmp = i.strip().split('\t')
		uniref = tmp[0]
		lineage = tmp[3].upper()
		uniref2lineage[uniref] = lineage

for i in open(args.path+'/marker_gene_group_lens.txt'):
	tmp = i.strip().split('\t')
	MG_group = tmp[0]
	MG_MSA_len = float(tmp[1])
	busco2protlen[MG_group] = MG_MSA_len

print "Extracted taxonomy information!!!"

# dictionary of marker genes reads aligned to
classifiers = {i.strip().split('\t')[1]:{} for i in open(args.d)}

print "Number of candidate marker genes found in sample: ",len(classifiers)

# collect read file information --> average read length and number of reads
for h,i in enumerate(open('read_file_info.txt')):
	if h == 0: read_count = float(i.strip())
	elif h == 1: av_read_len = float(i.strip())

# collect uniref to marker gene group mappings and the structure of gaps for each 
# uniref per their multiple sequence alignment with marker gene groups
marker_gene_groups = {}
marker_gene_indices = {}

for i in open(args.path+'/global_msa_indices.txt'):
	tmp = i.strip().split('\t')
	uniref = tmp[0]
	if uniref not in classifiers: continue
	F_name = tmp[1]
	F_gap_starts = map(int,tmp[2].split(','))
	F_gap_lens = map(int,tmp[3].split(','))
	G_name = tmp[4]
	G_gap_starts = map(int,tmp[5].split(','))
	G_gap_lens = map(int,tmp[6].split(','))
	S_name = tmp[7]
	S_gap_starts = map(int,tmp[8].split(','))
	S_gap_lens = map(int,tmp[9].split(','))
	marker_gene_groups[F_name] = 0
	marker_gene_groups[G_name] = 0
	marker_gene_groups[S_name] = 0
	marker_gene_indices[uniref] = {F_name:[F_gap_starts,F_gap_lens],G_name:[G_gap_starts,G_gap_lens],S_name:[S_gap_starts,S_gap_lens]}

print "Extracted MSA groups and coordinates!!!"

# collect coefficients of spline and logistic functions, as well as classification power by region, per marker gene group
# format of slopes file
# Marker_gene_group     Knots   A       B       C       D       LR_slopes       LR_y_ints       Classification_power
for h,i in enumerate(open(args.path+'/slopes2.threshold.logistic.txt')):
	if h == 0: continue
	tmp = i.strip().split('\t')
	marker_gene_group = tmp[0]
	knots = map(float,tmp[1].split(','))
	A = map(float,tmp[2].split(','))
	B = map(float,tmp[3].split(','))
	C = map(float,tmp[4].split(','))
	D = map(float,tmp[5].split(','))
	LR_slopes = [float(i) if i != 'NA' else 'NA' for i in tmp[6].split(',')]
	LR_y_ints = [float(i) if i != 'NA' else 'NA' for i in tmp[7].split(',')]
	Classification_power = map(float,tmp[8].split(','))
	if marker_gene_group in marker_gene_groups:
		marker_gene_groups[marker_gene_group] = [knots,A,B,C,D,LR_slopes,LR_y_ints,Classification_power]

print "Extracted marker gene classifiers!!!"

# classify each read per marker gene group it aligns to
with open('reads_classified.txt','w') as out:
	for i in open(args.d):
		tmp = i.strip().split('\t')
		query = tmp[0]
		uniref = tmp[1]
		aln_len = int(tmp[3])
		qstart = int(tmp[8])
		bitscore = float(tmp[11])
		qslope = bitscore/aln_len
		if uniref not in marker_gene_indices: continue
		results = {}
		for marker_gene_group,gaps in marker_gene_indices[uniref].items():
			adjusted_start = adjust_start(qstart,gaps)
			if marker_gene_groups[marker_gene_group] == 0: continue
			for knot,A,B,C,D,LR_slope,LR_y_int,Classification_power in zip(*marker_gene_groups[marker_gene_group]):
				dist_from_knot = abs(adjusted_start - knot)
				if dist_from_knot <= 10:
					probability,classification_power_of_region = classify(dist_from_knot,qslope,A,B,C,D,LR_slope,LR_y_int,Classification_power)
					taxa_rank = nodes[marker_gene_group.split('_')[1]]
					results[taxa_rank]= [marker_gene_group,probability,classification_power_of_region]
					break
		if results == {}: continue
		for rank in ['species','genus','family']:
			if rank in results:
				out.write(query+'\t'+uniref+'\t'+results[rank][0]+'\t'+str(results[rank][1])+'\t'+str(results[rank][2])+'\n')
				break


print "Finished classifying individual reads!!! Aggregating results..."

taxa_counts = {}
best_probs = {}

with open('Taxonomic_report.txt','w') as out:
	for i in open('reads_classified.txt'):
		tmp = i.strip().split('\t')
		query = tmp[0]
		uniref = tmp[1]
		taxaID = tmp[2].split('_')[1]
		lineage_trimmed = fullLineage[taxaID]
		r_prob = float(tmp[3])
		classification_power = float(tmp[4])
		r_prob = r_prob*classification_power
		if r_prob >= (float(args.th)):
			if query not in best_probs:
				best_probs[query] = [lineage_trimmed]
			else:
				best_probs[query].append(lineage_trimmed)
	for query,lines in best_probs.items():
		lines = list(set(lines))
		if len(lines) > 1:
			lca = LCA(lines)
		else:
			lca = map(lambda x:x.strip(),lines[0].split(';'))
		if lca == []: continue
		line = ''
		for i in lca:
			if i == '': continue
			line += i+';'
			if line not in taxa_counts:
				taxa_counts[line] = 1
			else:
				taxa_counts[line] += 1
	for k,v in sorted(taxa_counts.items()):
		relative_abundance = str((v*av_read_len)/(busco2protlen[names[k.split(';')[-2]]]*read_count))
		out.write(k+'\t'+str(v)+'\t'+relative_abundance+'\n')


print "Taxonomic report has successfully been written!!! Analysis completed successfully!"



