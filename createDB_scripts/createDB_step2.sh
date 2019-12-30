#!/bin/sh 
#$ -N step2
#$ -j y 
#$ -pe mpi 1
#$ -cwd

# To run: qsub -tc 50 -t 1-947 createDB_step2.sh ../data/fastas.txt


# Load modules and set paths
#
PWD=`pwd`
diamond="/nfs/sw/apps/diamond/diamond"
kaiju="/nfs/sw/apps/kaiju/"
createDB=$PWD
data=$PWD"/../data/"
busco=$PWD"/../busco/"
run_pipeline=$PWD"/../run_pipeline_scripts/"


cd $data


# Create array 
#
INPUTFILE=$1 
INPUTS=$(awk "NR==$SGE_TASK_ID" $INPUTFILE) 
INPUTSarray=($INPUTS)


# Process input 
#
FILE1=${INPUTSarray[0]}


# Align k-mers to eukaryotic busco proteins
#
# $diamond blastp --min-score 50 --db euk_busco --query $FILE1 --threads 10 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq_gapped sseq_gapped --out $FILE1.diamond --unal 0


# Find taxonomic labels, alignment slopes and start position of aligned k-mers; then sort files
#
python $createDB/compare_alns.py -d $FILE1.diamond -f euk_functions.txt -i global_msa_indices.txt -names names -nodes nodes -o $FILE1.diamond.slopes2
sort $FILE1.diamond.slopes2 > $FILE1.diamond.slopes2.sort



