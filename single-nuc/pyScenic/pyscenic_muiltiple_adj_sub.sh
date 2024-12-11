#!/bin/sh

# Where PBS script output and error files are stored.
PBS_OUTPUT="$RDS/home/pyScenic"

mkdir -p $PBS_OUTPUT
cd $PBS_OUTPUT


# Path to the .pbs script to run arboreto
PBS_SCRIPT="sub_pyScenic_adj_multiple.sh"

END=10
for i in $(seq 1 $END)
do
    dataset="global_all_stringent_sub$i"
    input=$dataset".loom"
    output=$dataset"adj.tsv"
	A=$(qsub -v INPUT=$input,OUTPUT=$output $PBS_SCRIPT) 
	echo $dataset $A 
done