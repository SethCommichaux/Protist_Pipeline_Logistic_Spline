#! /bin/sh
#$ -N step3
#$ -j y
#$ -pe mpi 1
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


# Merge sort diamond slopes files
#
# sort -m -s -T $data *slopes2.sort > slopes2.txt.sorted


# fit decision threshold and logistic functions to each marker gene
#
python $createDB/minimize_error_logistic.py -s slopes2.txt.sorted -o slopes2.threshold.logistic.txt

