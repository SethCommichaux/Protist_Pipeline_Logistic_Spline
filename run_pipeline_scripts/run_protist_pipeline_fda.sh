#!/bin/sh
#$ -N proteinDB
#$ -j y
#$ -pe mpi 12
#$ -cwd


# Load modules and software paths into environment
#
module load java
module load biopython
module load pandas
module load trinity
diamond="/nfs/sw/apps/diamond/diamond"
kaiju="/nfs/sw/apps/kaiju/"
createDB="/lustre/projects/SethCommichaux/Busco_Protist_Pipeline/createDB_scripts/"
protist_data="/lustre/projects/SethCommichaux/Busco_Protist_Pipeline/data/"
run_pipeline="/lustre/projects/SethCommichaux/Busco_Protist_Pipeline/run_pipeline_scripts/"
kaijuDB="/lustre/projects/SethCommichaux/Busco_Protist_Pipeline/data/euk_busco.pep.kaiju.fmi"
queryDB="/lustre/projects/SethCommichaux/Busco_Protist_Pipeline/data/euk_busco"
trimmomatic="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/Trimmomatic-0.38/trimmomatic-0.38.jar"
adapters="/lustre/projects/SethCommichaux/ProteinDB_protist_pipeline/Trimmomatic-0.38/adapters/TruSeq3-SE.fa"


# Input and output directories
#
input=`pwd`
output=`pwd`
mkdir $output


################################################################
################################################################
################################################################

for i in $input/*fastq

do time

# Fastq file(s) to be analyzed
#
#read_count=`grep -c "^+" $i`
#reads_fastq=$(basename ${i%.fastq}.trimmed.fastq)
reads_fastq=$i

# Output file
#
out=$output/$(basename ${i%.fastq})
mkdir $out
cd $out


# Trim/Filter raw reads
#
# java -jar $trimmomatic SE -threads 12 $i $reads_fastq ILLUMINACLIP:$adapters:2:30:10 MAXINFO:120:0.2


# Run kaiju to query fastq reads against protein sequence binning databse (binningDB.fasta)
#
$kaiju/kaijux -f $kaijuDB -i $in/$reads_fastq -z 12 -m 9 | grep "^C" > kaiju


# Extract reads that aligned to binning database
#
python $run_pipeline/extract_kaiju_reads.py -k kaiju -s $in/$reads_fastq -o kaiju.fasta


# Align binned reads, with Diamond, to queryDB
#
# --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen
time $diamond blastx --sensitive --evalue 0.0000000001 --db $queryDB --query kaiju.fasta --threads 12 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen --out kaiju.fasta.diamond


# Process diamond output
#
python $run_pipeline/subset_diamond_best_bitscore.py -d kaiju.fasta.diamond -o kaiju.fasta.diamond.subset
python $run_pipeline/classify_spline_logistic.py -d kaiju.fasta.diamond.subset -path $protist_data -th 0.9

##############################################################################################################
# Assemble binned reads with Trinity
#
# Trinity --seqType fa --single kaiju.fasta --max_memory 100G --output trinity --CPU 12


# Classify assembled reads
#
# time $diamond blastx --sensitive --evalue 0.0000000001 --db $queryDB --query trinity/Trinity.fasta --threads 12 --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq qlen --out Trinity.fasta.diamond


# Process diamond output for assembly
#
# python $run_pipeline/subset_diamond_best_bitscore.py -d Trinity.fasta.diamond -o Trinity.fasta.diamond.subset
# python $run_pipeline/score_logistic.py -d Trinity.fasta.diamond.subset -s $protist_data/sparsified_test.txt -o Taxonomic_report_trinity.txt -f $protist_data/euk_functions.txt -no $protist_data/nodes -na $protist_data/names
##############################################################################################################


# Create files for producing sankey diagrams
#
# python $run_pipeline/sankey_preprocess.py -nodes $protist_data/nodes.dmp -lineage $protist_data/fullnamelineage.dmp -results $out/kaiju.fasta.diamond.subset.results -o $out


done


# Create merged abundance output tables
#
# cd $output
# python $run_pipeline/merge_tables.py $output
# python $run_pipeline/filter_merged_abundance_table_protist.py -nodes $protist_data/nodes -names $protist_data/names -t merged_abundance_table.txt



