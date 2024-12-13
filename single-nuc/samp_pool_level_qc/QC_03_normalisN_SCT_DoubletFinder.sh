################################################################################################

#!/bin/bash
#PBS -l select=1:ncpus=12:mem=100gb:ngpus=4:gpu_type=RTX6000
#PBS -l walltime=6:00:00
#PBS -N seurat_norm_doub_it4

tid=$PBS_ARRAY_INDEX
module load anaconda3/personal
source activate snRNAseq

time R --vanilla --args ${tid} << "EOF"
args=commandArgs(trailingOnly = TRUE)
tid=as.numeric(as.character(args[1]))

library(SingleCellExperiment) 
library(Seurat) 
library(DoubletFinder)

## Step 1: Seurat normalisation
## Step 2: Run DoubletFinder

# load files
cleanFiles = list.files(.../01_quality_filter/",full.names = T)

# load genotype annotations per cell per sample
geno_donor_id = readRDS(".../donor_geno_anno.RDS")

# array job
njobs=18
interval=ceiling(length(cleanFiles)/njobs)
start=(interval*(tid-1))+1
end=ifelse((interval*tid)<length(cleanFiles), (interval*tid), length(cleanFiles))
cleanFiles=cleanFiles[start:end]

# loop across samples
for(c in 1:length(cleanFiles)){
  
  print(paste("normalising sample",basename(cleanFiles[c]),sep=" "))
  clean = readRDS(cleanFiles[c])
  print(dim(clean))

  # add genotype doublet metadata

  sample = unlist(strsplit(basename(cleanFiles[c]),split="_lib_clean.rds"))
  donor_gt = geno_donor_id[[sample]]
  metadata = clean@meta.data
  metadata$donor_id = donor_gt$donor_id[match(metadata$cells,donor_gt$cell)]
  metadata$GT = ifelse(metadata$donor_id=="doublet","Doublet","Singlet")
  clean@meta.data = metadata 
  
  doublet_fraction = sum(metadata$donor_id == "doublet") / length(metadata$donor_id)
  
  # normalise

  clean = NormalizeData(clean, verbose = FALSE)
  
  # scale data (old seurat method, as may use for visualisation)
  
  clean <- ScaleData(clean, features = rownames(clean), vars.to.regress = "percent.mt" )
  
  # SCT v2 with regularisation
  
  clean = SCTransform(clean,  vst.flavor = "v2", vars.to.regress = "percent.mt")
  #clean = SCTransform(clean, vars.to.regress = "percent.mt")
  
  saveRDS(clean,file.path("/rds/general/project/lms-scott-raw/live/snRNA_seq_processing/human_WAT/Phase_4/seurat/it4/02_SCT_transform",basename(cleanFiles[c])))
  
  clean = RunPCA(clean, verbose = FALSE)
  clean = RunUMAP(clean, dims = 1:30, verbose = FALSE)
  clean = RunTSNE(clean, dims = 1:30, verbose = FALSE)
  clean = FindNeighbors(clean, dims = 1:30, verbose = FALSE)
  clean = FindClusters(clean, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0), verbose = FALSE)

  saveRDS(clean,file.path("/rds/general/project/lms-scott-raw/live/snRNA_seq_processing/human_WAT/Phase_4/seurat/it4/02_SCT_transform",basename(cleanFiles[c])))
  
  sn_lib = clean
  rm(clean)
  
  set.seed(2020)
  
  ## pK Identification (with ground-truth) ---------------------------------------------------------------------------------------
  
  sweep.res.doublet <- suppressWarnings( paramSweep_v3(sn_lib, PCs = 1:30, sct = T) )

  gt.calls <- sn_lib@meta.data[rownames(sweep.res.doublet[[1]]), "GT"]

  sweep.stats.doublet <- summarizeSweep(sweep.res.doublet, GT = TRUE, GT.calls = gt.calls )

  bcmvn.doublet <- find.pK(sweep.stats.doublet)
  
  tiff(paste0("/rds/general/project/lms-scott-raw/live/snRNA_seq_processing/human_WAT/Phase_4/seurat/it4/03_DoubletFinder/",basename(cleanFiles[c]),".tiff"))
  bcmvn.doublet
  dev.off()
  
  pK <- bcmvn.doublet$pK[which.max(bcmvn.doublet$MeanAUC)]
  pK <- as.numeric(levels(pK))[pK]
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  annotations <- sn_lib@meta.data$SCT_snn_res.0.8
  homotypic.prop <- modelHomotypic(annotations)
  
  nExp_0.05 <- round(0.05*nrow(sn_lib@meta.data))
  nExp_0.05_adj <- round(nExp_0.05*(1-homotypic.prop))
  nExp_0.075 <- round(0.075*nrow(sn_lib@meta.data))
  nExp_0.075_adj <- round(nExp_0.075*(1-homotypic.prop))
  nExp_0.1 <- round(0.1*nrow(sn_lib@meta.data))  
  nExp_0.1_adj <- round(nExp_0.1*(1-homotypic.prop))
  nExp_0.2 <- round(0.2*nrow(sn_lib@meta.data)) 
  nExp_0.2_adj <- round(nExp_0.2*(1-homotypic.prop))
  
  # Run using genotype doublet estimates 

  nExp_GT <- round(doublet_fraction*nrow(sn_lib@meta.data)) 
  nExp_GT_adj <- round(nExp_GT*(1-homotypic.prop))

  sn_lib <-  suppressWarnings(doubletFinder_v3(sn_lib, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_0.05, reuse.pANN = FALSE, sct = T))
  reuse.pANN = paste("pANN",0.25,pK,nExp_0.05,sep="_")
  
  sn_lib <- doubletFinder_v3(sn_lib, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_0.075, reuse.pANN = reuse.pANN, sct = T)
  sn_lib <- doubletFinder_v3(sn_lib, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_0.1, reuse.pANN = reuse.pANN, sct = T)
  sn_lib <- doubletFinder_v3(sn_lib, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_0.2, reuse.pANN = reuse.pANN, sct = T)
  sn_lib <- doubletFinder_v3(sn_lib, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_GT, reuse.pANN = reuse.pANN, sct = T)
  
  sn_lib <- doubletFinder_v3(sn_lib, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_0.05_adj, reuse.pANN = reuse.pANN, sct = T)
  sn_lib <- doubletFinder_v3(sn_lib, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_0.075_adj, reuse.pANN = reuse.pANN, sct = T)
  sn_lib <- doubletFinder_v3(sn_lib, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_0.1_adj, reuse.pANN = reuse.pANN, sct = T)
  sn_lib <- doubletFinder_v3(sn_lib, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_0.2_adj, reuse.pANN = reuse.pANN, sct = T)
  sn_lib <- doubletFinder_v3(sn_lib, PCs = 1:10, pN = 0.25, pK = pK, nExp = nExp_GT_adj, reuse.pANN = reuse.pANN, sct = T)
  
  saveRDS(.../02_DoubletFinder/new",basename(cleanFiles[c])))
  
}

"EOF"
