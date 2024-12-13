
#!/bin/sh
#PBS -S /bin/bash
## PBS -N v.$spl
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -j oe
# This version works with or withour genotype. Also takes --callAmbiantRNA
# Some large sample need more mem and walltime at noGT mode.

source $HOME/.bashrc
module load anaconda3/personal
conda activate base

if [ -z ${spl+x} ]; then
    spl=$1
    basedir=$2
    n_donor=$3
fi
if [ -z ${indir+x} ]; then	# this is phase 4 default
	indir=$basedir/cellSNP-lite/G38contig_order/$spl
fi
if [ -z ${outdir+x} ]; then	# this is for ambiant RNA call with impGT
	outdir=$basedir/vireo/with_geno/impGT.aRNA/$spl
fi
if [ -z ${donor_GT+x} ]; then
	donor_op=""
else
	donor_op="-d $donor_GT --genoTag GT"
	ls -la $donor_GT
fi
if [ -z ${n_donor+x} ]; then
	n_donor=4
fi
if [ -z ${aRNA+x} ]; then
	aRNA_op=""
else
	aRNA_op="--callAmbientRNAs"
fi
if [ -z ${x_donor+x} ]; then
	x_donor_op=""
else
	x_donor_op="--extraDonor $x_donor"
fi
lsr=/rds/general/project/lms-scott-raw/live
export PATH=$lsr/resources/software/vireo/bin:$lsr/resources/software/samtools/bcftools-1.10.2/bin:$lsr/resources/software/samtools/htslib-1.10.2:$PATH
export LD_LIBRARY_PATH=$lsr/resources/software/samtools/htslib-1.10.2:$LD_LIBRARY_PATH
export PYTHONPATH=$lsr/resources/software/vireo/:$PYTHONPATH

cd $basedir
echo -n 'starting at: '; date   # debug
#CELL_DATA=$basedir/cellSNP-lite/$spl
CELL_DATA=$indir
echo "input dir: $CELL_DATA"
ls -la $CELL_DATA # $outdir

vireo -c $CELL_DATA $donor_op -N $n_donor -o $outdir -p 2 $aRNA_op $x_donor_op
