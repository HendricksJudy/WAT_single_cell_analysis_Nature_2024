
#!/bin/sh
#PBS -S /bin/bash
## PBS -N CSl
#PBS -l walltime=8:00:00
#PBS -l select=1:ncpus=16:mem=24gb
#PBS -j oe

source $HOME/.bashrc
module load anaconda3/personal
conda activate vireo

lsr=/rds/general/project/lms-scott-raw/live
export PATH=$lsr/resources/software/samtools/samtools-1.10/bin:$lsr/resources/software/samtools/bcftools-1.10.2/bin:$lsr/resources/software/samtools/htslib-1.10.2:$PATH
export LD_LIBRARY_PATH=$lsr/resources/software/samtools/htslib-1.10.2

if [ -z ${spl+x} ]; then
    spl=$1
    inbam=$2	# $lsr/combined_ngs_cellranger_outs_human_WAT_obese_lean_wl/Processed_combined_introns/
    #output_base=$3	# $lsr/snRNA_seq_processing/human_WAT/Phase_4/
    #REGION_VCF=$4
    #outdir=$5
    #MAFcut=$6
fi
spl_=`echo $spl|sed -e 's/-//'`	# compatibility for Phase 3 sample name
if [ -z ${inbam+x} ]; then	# Phase 4 default except WS1 (Pool1B)
    inbam=$lsr/combined_ngs_cellranger_outs_human_WAT_obese_lean_wl/Processed_combined_introns/$spl/outs/possorted_genome_bam.bam
fi
if [ -z ${output_base+x} ]; then
	output_base=$lsr/snRNA_seq_processing/human_WAT/Phase_4
fi
if [ -z ${outdir+x} ]; then
	outdir=$output_base/cellSNP-lite/cb_bcode/$spl_
fi
if [ -z ${REGION_VCF+x} ]; then
	#REGION_VCF=./resources/genome1k_from_vireo/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.srt.vcf.gz
    REGION_VCF=$lsr/resources/genome1k_from_vireo/genome1K.phase3.SNP_AF5e4.chr1toX.GRCh38.srt.vcf.gz
fi
if [ -z ${MAFcut+x} ]; then
	MAFcut=0.05 # 0.1 is default but we use 0.05
fi
if [ -z ${bcode+x} ]; then
    bcode=$output_base/intermediate_input/$spl_/cb_barcodes.csv
fi

cd $output_base
echo -n 'starting at: '; date   # debug
ls -la $inbam $bcode $REGION_VCF # $outdir

cellsnp-lite -s $inbam -b $bcode -O $outdir -R $REGION_VCF -p 16 --minMAF $MAFcut --minCOUNT 20 --gzip
