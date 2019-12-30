import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-bf", help="Path to busco fungi ancestral file")
parser.add_argument("-be", help="Path to busco eukaryota ancestral file")
parser.add_argument("-bp", help="Path to busco protists ancestral file")
parser.add_argument("-ep", help="Path to busco search eukaryota to protists file")
parser.add_argument("-ef", help="Path to busco search eukaryota to fungi file")
parser.add_argument("-o", help="Desired name for output file")
args = parser.parse_args()

# example command
# python map_buscos.py -bf fungi_odb9/ancestral -be eukaryota_odb9/ancestral -bp protists_ensembl/ancestral -ep run_protist_ancestral/full_table_protist_ancestral.tsv -ef run_fungi_ancestral/full_table_fungi_ancestral.tsv -o busco2busco_mappings.txt

# retrieve fasta identifiers of ancestral, busco genes 
fungi = {i.strip().replace(">",""):0 for i in open(args.bf) if i[0] == ">"}
eukaryota = {i.strip().replace(">",""):[] for i in open(args.be) if i[0] == ">"}
protists = {i.strip().replace(">",""):0 for i in open(args.bp) if i[0] == ">"}

# retrieve busco hmm results mapping eukaryota ancestral variants to protist and fungi variants
euk2protists = {i.strip().split('\t')[0]:i.strip().split('\t')[2] for i in open(args.ep) if i[0] != "#" if i.strip().split('\t')[1] == "Complete"}
euk2fungi = {i.strip().split('\t')[0]:i.strip().split('\t')[2] for i in open(args.ef) if i[0] != "#" if i.strip().split('\t')[1] == "Complete"}

# collect mappings between eukaryota ancestral variants to protist and fungi variants
for k,v in eukaryota.items():
	if k in euk2fungi:
		eukaryota[k].append(euk2fungi[k])
	if k in euk2protists:
		eukaryota[k].append(euk2protists[k])


r = {}

with open(args.o,'w') as out:
	for k in fungi: # output fungi marker gene identifiers with no eukaryota marker gene equivalents
		if k not in euk2fungi.values():	
			out.write(k+'\n')
	for k in protists: # output protist marker gene identifiers with no eukaryota marker gene equivalents
		if k not in euk2protists.values():
			out.write(k+'\n')
	for k,v in eukaryota.items(): # output all eukaryota identifiers, including their mappings to protist/fungi equivalents
		if v == []:
			out.write(k+'\n')
		else:
			out.write(k+'\t'+'\t'.join(v)+'\n')
	





