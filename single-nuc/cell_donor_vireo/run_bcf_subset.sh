
#!/bin/sh
#PBS -S /bin/bash
## PBS -N subset_$spl
#PBS -l walltime=4:00:00
#PBS -l select=1:ncpus=8:mem=16gb
#PBS -j oe

export PATH=/rds/general/project/lms-scott-raw/live/resources/software/samtools/bcftools-1.10.2/bin:/rds/general/project/lms-scott-raw/live/resources/software/samtools/htslib-1.10.2:$PATH
export LD_LIBRARY_PATH=/rds/general/project/lms-scott-raw/live/resources/software/samtools/htslib-1.10.2

lsr=/rds/general/project/lms-scott-raw/live
if [ -z ${spl+x} ]; then
    spl=$1
    REGION_VCF=$2
    GT_VCF=$3
fi
if [ -z ${outdir+x} ]; then
	outdir=$lsr/snRNA_seq_processing/human_WAT/Phase_3/vireo/wtGT/imputed/$spl
fi
if [ -z ${REGION_VCF+x} ]; then
	REGION_VCF=$lsr/resources/genome1k_from_vireo/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.srt.vcf.gz
fi
if [ -z ${GT_VCF+x} ]; then
	GT_VCF=$lsr/Adipocyte_GWAS/Geno_combined/imputed_genotype/G38/imputed_G38_snRNAseq_gt.vcf.gz
fi
echo -n "Starting: ";date
ls -la $outdir $REGION_VCF $GT_VCF
cd $lsr
bcftools view $GT_VCF -S $outdir/donor.subspl.lis -R $REGION_VCF -v snps -I -Oz -o $outdir/donor.sub.vcf.gz --threads 8
echo -n "$spl subset done: "; date

bcftools index --threads 8 $outdir/donor.sub.vcf.gz
echo -n "idx done: "; date
ls -lat $outdir
