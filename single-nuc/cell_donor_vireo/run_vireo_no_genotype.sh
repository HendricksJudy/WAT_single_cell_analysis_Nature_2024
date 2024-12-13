#!/bin/sh
#PBS -S /bin/bash
## PBS -N vireoNoGT_$spl
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=2:mem=8gb
#PBS -j oe

source $HOME/.bashrc
module load anaconda3/personal
conda activate base

if [ -z ${spl+x} ]; then
    spl=$1
    basedir=$2
    n_donor=$3
fi
if [ -z ${outdir+x} ]; then
	outdir=$basedir/vireo/noGT/$spl
fi
if [ -z ${n_donor+x} ]; then
	n_donor=4
fi
if [ -z ${x_donor+x} ]; then
	x_donor_op=""
else
	x_donor_op="--extraDonor $x_donor"
fi
lsr=/rds/general/project/lms-scott-raw/live
export PATH=$lsr/resources/software/vireo/bin:$lsr/resources/software/samtools/bcftools-1.10.2/bin:$lsr/resources/software/samtools/htslib-1.10.2:$PATH
export LD_LIBRARY_PATH=$lsr/resources/software/samtools/htslib-1.10.2
export PYTHONPATH=$lsr/resources/software/vireo/

cd $lsr
echo -n 'starting at: '; date   # debug
CELL_DATA=$basedir/cellSNP-lite/$spl
#DONOR_GT_F=./snRNA_seq_processing/human_WAT/Phase_2/vireo/donor_GT/donor_GT.vcf.gz
ls -la $CELL_DATA # $outdir $DONOR_GT_F

vireo -c $CELL_DATA -N $n_donor -o $outdir -p 2 $x_donor_op
