################################################################################################

#!/bin/bash
#PBS -l select=1:ncpus=48:mem=120gb
#PBS -l walltime=24:00:00
#PBS -N Nebula_MO_NW_rsn_cont

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

#load seurat object all samples 
swat = readRDS(".../swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.rds")

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

anno = "cell_type_am_fine"
#anno = "cell_state_am_long"

cell_types = names(table(swat@meta.data[[anno]]))
cell_types 

# set categorical covariates as factors

coldat = swat@meta.data
coldat$eth = as.factor(coldat$eth)
coldat$sex = as.factor(coldat$sex)
coldat$dataset = as.factor(coldat$dataset)

swat@meta.data = coldat

# set directory

dir = "/.../nebula/level1/"
#dir = "/.../nebula/level2/"

# run nebula

njobs=13
interval=ceiling(length(cell_types)/njobs)
start=(interval*(tid-1))+1
end=ifelse((interval*tid)<length(cell_types), (interval*tid), length(cell_types))
cell_types=cell_types[start:end]

for(c in 1:length(cell_types)){
  
  print(cell_types[c])
  
  res = nebula_diff_gex(data = swat, idents = anno, celltype = cell_types[c], comparison = c("Lean","Obese"), covariates = c("mt.percent","ribo.percent","age","sex","eth"), offset = "sizeFactors", sample_id ="sample", downsample = FALSE, method = "LN")
    
  res2 = cbind(res$summary,res$overdispersion, res$convergence, res$algorithm)
  res2 = merge(res2, res$pct, by = "gene", all.x=T)
  res2 = res2[order(res2$p_conditionObese),]

  #colnames(res2) = c(colnames(res2)[1:11],"subject_overdisp","cell_overdisp","convergence","algorithm",colnames(res2)[16:18])
  print(head(res2))

  fname = paste("obese",anno,cell_types[c],"LN","nebula_diff_expn_tech_bio_low_mt_rb_cov.RDS",sep="_")
  saveRDS(res,file = file.path(dir,fname))
  
  fname2 = paste("obese",anno,cell_types[c],"LN","nebula_diff_expn_tech_bio_low_mt_rb_cov.csv",sep="_")
  write.csv(res2,file = file.path(dir,fname2))
  
}

'EOF'

