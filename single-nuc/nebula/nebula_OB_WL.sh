################################################################################################

#!/bin/bash
#PBS -l select=1:ncpus=48:mem=124gb
#PBS -l walltime=24:00:00
#PBS -N Nebula_MO_WL_pair

tid=$PBS_ARRAY_INDEX
module load anaconda3/personal
source activate nebula

time R --vanilla --args ${tid} << "EOF"
args=commandArgs(trailingOnly = TRUE)
tid=as.numeric(as.character(args[1]))

#################################################################################################

library(Seurat)
source(".../Nebula_script.sh")

#### Nebula test run

# Run in cluster

swat = readRDS("/.../nebula/swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.rds")

# remove ribo and mt <1% 

swat = subset(swat, subset = mt.percent < 1)
swat = subset(swat, subset = ribo.percent < 1)

# remove emont samples

swat = subset(swat,  subset = dataset == "scott")

# calcuate size scaling factors

sce = Seurat::as.SingleCellExperiment(swat, assay = "RNA")
sizeFactors = scater::librarySizeFactors(sce)
swat$sizeFactors = sizeFactors

colnames(swat@meta.data)

# set annotation level

#anno = "cell_type_am_fine"
anno = "cell_state_am_long"

cell_types = names(table(swat@meta.data[[anno]]))
cell_types 

# set categorical covariates as factors

coldat = swat@meta.data
coldat$eth = as.factor(coldat$eth)
coldat$sex = as.factor(coldat$sex)

swat@meta.data = coldat

# run nebula

#dir = "/rds/general/project/lms-scott-raw/live/snRNA_seq_processing/human_WAT/Phase_5/nebula/level1/"
dir = "/rds/general/project/lms-scott-raw/live/snRNA_seq_processing/human_WAT/Phase_5/nebula/level2/"

njobs=13
interval=ceiling(length(cell_types)/njobs)
start=(interval*(tid-1))+1
end=ifelse((interval*tid)<length(cell_types), (interval*tid), length(cell_types))
cell_types=cell_types[start:end]

for(c in 1:length(cell_types)){
  
  print(cell_types[c])
  
  res = nebula_diff_gex(data =  swat, idents = anno, celltype = cell_types[c], comparison = c("Obese","Weightloss"), covariates = c("mt.percent","ribo.percent"), offset = "sizeFactors", sample_id ="sample", downsample = FALSE, method = "LN")
  
  res2 = cbind(res$summary,res$overdispersion, res$convergence, res$algorithm)
  res2 = merge(res2, res$pct, by = "gene", all.x=T)
  res2 = res2[order(res2$p_conditionWeightloss),]
  #colnames(res2) = c(colnames(res2)[1:11],"subject_overdisp","cell_overdisp","convergence","algorithm",colnames(res2)[16:18])
  print(head(res2))

  fname = paste("weightloss",anno,cell_types[c],"LN","nebula_diff_expn_tech_paired_low_mt_rb.RDS",sep="_")
  saveRDS(res,file = file.path(dir,fname))
  
  fname2 = paste("weightloss",anno,cell_types[c],"LN","nebula_diff_expn_tech_paired_low_mt_rb.csv",sep="_")
  write.csv(res2,file = file.path(dir,fname2))
  
}
 
'EOF'

