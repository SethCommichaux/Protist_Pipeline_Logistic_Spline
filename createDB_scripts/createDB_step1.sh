#! /bin/sh
#$ -N DB_create
#$ -j y
#$ -pe mpi 12
#$ -cwd


# Load modules and software paths into environment
#
module load biopython
module load pandas
module load mafft
PWD=`pwd`
diamond="/nfs/sw/apps/diamond/diamond"
kaiju="/nfs/sw/apps/kaiju/"
createDB=$PWD
data=$PWD"/../data/"
busco=$PWD"/../busco/"
run_pipeline=$PWD"/../run_pipeline_scripts/"


# Change to data directory
#
mkdir $data
cd $data


# Download, decompress NCBI taxonomy files
#
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
# tar xvf new_taxdump.tar.gz
# grep -i 'scientific name' names.dmp | awk -F "\t|\t" '{print $1 "\t" $3}' > names
# awk -F "\t|\t" '{print $1 "\t" $5}' nodes.dmp > nodes
# rm nodes.dmp names.dmp new_taxdump.tar.gz rankedlineage.dmp taxidlineage.dmp type* citations.dmp delnodes.dmp division.dmp gencode.dmp host.dmp merged.dmp


# Download uniprot idmapping file and UniRef100 protein sequence file
#
# wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
# wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
# gunzip *gz


# Download BUSCO datasets for eukaryota, fungi and protists
#
# wget http://busco.ezlab.org/v2/datasets/eukaryota_odb9.tar.gz
# wget http://busco.ezlab.org/v2/datasets/fungi_odb9.tar.gz
# wget http://busco.ezlab.org/v2/datasets/protists_ensembl.tar.gz
# tar xvf eukaryota_odb9.tar.gz
# tar xvf fungi_odb9.tar.gz
# tar xvf protists_ensembl.tar.gz
# rm *.tar.gz


# Identify homologs between BUSCO eukaryota, fungi and protist datasets
#
# python $busco/scripts/run_BUSCO.py -i protists_ensembl/ancestral -o protist_ancestral -c 12 -m prot -l eukaryota_odb9/
# python $busco/scripts/run_BUSCO.py -i fungi_odb9/ancestral -o fungi_ancestral -c 12 -m prot -l eukaryota_odb9/
# python $createDB/map_buscos.py -bf fungi_odb9/ancestral -be eukaryota_odb9/ancestral -bp protists_ensembl/ancestral -ep run_protist_ancestral/full_table_protist_ancestral.tsv -ef run_fungi_ancestral/full_table_fungi_ancestral.tsv -o busco2busco_mappings.txt


# Extract eukarota and non-eukaryota protein sequences; identify protist and fungi marker genes
#
# python $createDB/splitUniRef100.py -g uniref100.fasta -t fullnamelineage.dmp
# python $busco/scripts/run_BUSCO.py -i fungi.pep -o eukaryota.busco.fungi -c 12 -m prot -l fungi_odb9/
# python $busco/scripts/run_BUSCO.py -i protist.pep -o eukaryota.busco.protist -c 12 -m prot -l protists_ensembl/
# cat fungi.pep protist.pep > eukaryota.pep
# python $busco/scripts/run_BUSCO.py -i eukaryota.pep -o eukaryota.busco.euk -c 12 -m prot -l eukaryota_odb9/
# grep -v "^#" */full_table* | cut -f 1,3 -d$'\t' | sort | uniq > busco_marker_gene_ids.txt
# python $createDB/extract_busco.py -b busco_marker_gene_ids.txt -e eukaryota.pep -o euk_busco.pep
# rm protist.pep fungi.pep eukaryota.pep


# Create annotation file for protist and fungi marker genes
#
# python $createDB/processDB.py -q euk_busco.pep -u idmapping.dat -t fullnamelineage.dmp -o euk_functions.txt
# rm idmapping.dat


# Split marker genes into taxonomic-based clusters
#
# mkdir taxa_groups
# cd taxa_groups
# mkdir singletons MSA
# python $createDB/split_into_groups.py -b2b ../busco2busco_mappings.txt -e ../euk_busco.pep -b ../busco_marker_gene_ids.txt -f ../euk_functions.txt -no ../nodes -na ../names
# cd MSA
# for i in *fasta; do mafft --thread 12 $i > $i.msa; done
# cd ../../


# Extract the global coordinates for each marker gene in each multiple sequence alignment
#
# python $createDB/global_msa_coordinates.py -nodes nodes -names names -o global_msa_indices.txt -p taxa_groups/MSA/
# python $createDB/global_msa_coordinates.py -nodes nodes -names names -o global_msa_indices.txt -p taxa_groups/singletons/


# make marker gene database for diamond and kaiju softwares
#
# $diamond makedb --in euk_busco.pep --db euk_busco --threads 12
# $kaiju/mkbwt -o euk_busco.pep.kaiju -n 12 -l 100000 euk_busco.pep
# $kaiju/mkfmi euk_busco.pep.kaiju
# rm euk_busco.pep.kaiju.bwt euk_busco.pep.kaiju.sa


# Identify candidates for non-protist/non-fungi negatives dataset
#
# $kaiju/kaijup -a greedy -e 3 -f euk_busco.pep.kaiju.fmi -i non_eukaryota.pep -z 12 -m 12 -s 50 | grep "^C" > noneuk2euk
# python $createDB/extract_kaiju_reads.py -k noneuk2euk -s non_eukaryota.pep -o negatives.pep


# Extract kmers from positive and negative datasets
#
# python $createDB/kmer_split.py -k 70 -s 10 -i negatives.pep -o negatives70_10.pep
# python $createDB/kmer_split.py -k 70 -s 10 -i euk_busco.pep -o euk_busco70_10.pep
# cat negatives70_10.pep euk_busco70_10.pep > all70_10.pep
# python $createDB/diamond_preprocess_kmers.py all70_10.pep
# ls *fasta > fastas.txt




