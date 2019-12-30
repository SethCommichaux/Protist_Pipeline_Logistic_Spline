import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("-q", help="Path to protists.pep protein fasta file")
parser.add_argument("-t", help="Path to NCBI fullnamelineage.dmp taxonomy file")
parser.add_argument("-u", help="Path to UniProt protein ID mapping file")
parser.add_argument("-o", help="Desired name for output file")
args = parser.parse_args()


# parse NCBI taxonomy fullnamelineage.dmp file for taxaIDs, taxaNames and taxonomic lineages
def build_taxa_dict(): 
	taxaID2_Name_Lineage = {}
	for i in open(args.t):
		tmp = i.strip().split('\t|\t')
		taxaID = tmp[0]
		taxaName = tmp[1].upper().strip('\t|\t')
		lineage = tmp[2].upper().strip('\t|\t')
		taxaID2_Name_Lineage[taxaID] = [taxaName,lineage+taxaName]
	print "Taxa dictionary built!!!"
	return taxaID2_Name_Lineage


def uniprot_annotations(ids):
	results = {}
	for i in open(args.u):
		tmp = i.strip().split('\t')
		uniprotID = tmp[0]
		category = tmp[1]
		categoryID = tmp[2]
		if uniprotID not in ids: continue
		if uniprotID not in results:
			results[uniprotID] = [0,0,0]
		elif category == 'KEGG':
			results[uniprotID][0] = categoryID
		elif category == 'KO':
			results[uniprotID][1] = categoryID
		elif category == 'eggNOG':
			results[uniprotID][2] = categoryID
	return results


# Build taxa dict
taxaID2_Name_Lineage = build_taxa_dict()


# retrieve queryDB uniprotIDs
ids = {i.split('_')[1].split(' ')[0] for i in open(args.q) if i.startswith(">")}


# retrieve annotations
results = uniprot_annotations(ids)


with open(args.o,'w') as out:
	out.write('UNIREF100_ID\tTAXANAME\tTAXAID\tLINEAGE\tPROTEIN_NAME\tPROTEIN_LENGTH\tKEGG\tKO\teggNOG\n')
	for i in SeqIO.parse(args.q,'fasta'):
		UniRef100 = str(i.id)
		uniprot = UniRef100.replace('UniRef100_','')
		d = str(i.description)
		prot_len = str(len(i.seq))
		taxaID = d.split('TaxID=')[1].split(' ')[0]
		if taxaID not in taxaID2_Name_Lineage: continue
		lineage = taxaID2_Name_Lineage[taxaID][1]
		taxaName = taxaID2_Name_Lineage[taxaID][0]
		proteinName = ' '.join(d.split(' ')[1:]).split('n=')[0]
		if uniprot not in results:
			KEGG,KO,eggNOG = '0','0','0'
		else:
			KEGG = str(results[uniprot][0])
			KO = str(results[uniprot][1])
			eggNOG = str(results[uniprot][2])
		out.write(UniRef100+'\t'+taxaName+'\t'+taxaID+'\t'+lineage+'\t'+proteinName+'\t'+prot_len+'\t'+KEGG+'\t'+KO+'\t'+eggNOG+'\n')


