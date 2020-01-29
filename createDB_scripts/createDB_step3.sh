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


# Change to diamond results directory
#
cd $data/diamond_results/

# Merge sort diamond slopes files
#
sort -m -s -T $data/diamond_results/ *slopes2.sort > $data/slopes2.txt.sorted

# Change to data directory
#
cd $data/

# fit spline and logistic functions to each marker gene
#
python $createDB/minimize_error_logistic.py -s slopes2.txt.sorted -o slopes2.threshold.logistic.txt

