################################################################################################

#!/bin/bash
#PBS -l select=1:ncpus=16:mem=60gb
#PBS -l walltime=24:00:00
#PBS -N seurat_re_norm

tid=$PBS_ARRAY_INDEX
module load anaconda3/personal
source activate new

time R --vanilla --args ${tid} << "EOF"
args=commandArgs(trailingOnly = TRUE)
tid=as.numeric(as.character(args[1]))

library(Seurat) 

## Load Emont dataset, remove MT and RB >5%, split, and SCT transform

# load emont dataset

rosen_swat = readRDS(".../rosen_human_sat_no_svf.rds")
coldat = rosen_swat@meta.data
print(colnames(coldat))

# Add categorical metadata 

coldat$dataset = "rosen"
coldat = within(coldat, { age.range <- NA # need to initialize variable
  age.range[age == "nan"] <- NA 
  age.range[age < 20] <- "<20"  
  age.range[age >= 20 & age < 30] <- "20-30"
  age.range[age >= 30 & age < 40] <- "30-40"
  age.range[age >= 40 & age < 50] <- "40-50"
  age.range[age >= 50 & age < 60] <- "50-60"
  age.range[age >= 60 & age < 70] <- "60-70"
  age.range[age >= 70 & age < 80] <- "70-80"
  age.range[age >= 80] <- ">80"
} )
coldat$age.range = ifelse(coldat$age=="nan", NA, coldat$age.range)
coldat$age.range = factor(coldat$age.range, ordered = T, levels = c("<20","20-30","30-40","40-50","50-60","60-70","70-80",">80"))
coldat$sex = ifelse(coldat$sex == "Male","M","F")
coldat = within(coldat, { bmi.range <- NA # need to initialize variable
  bmi.range[bmi < 20] <- "<20"  
  bmi.range[bmi >= 20 & bmi < 25] <- "20-25"
  bmi.range[bmi >= 25 & bmi < 30] <- "25-30"
  bmi.range[bmi >= 30 & bmi < 35] <- "30-35"
  bmi.range[bmi >= 35 & bmi < 40] <- "35-40"
  bmi.range[bmi >= 40 & bmi < 45] <- "40-45"
  bmi.range[bmi >= 45 & bmi < 50] <- "45-50"
  bmi.range[bmi >= 50] <- ">50"
} )
coldat$bmi.range = factor(coldat$bmi.range, ordered = T, levels = c("20-25", "25-30", "30-35", "35-40", "40-45","45-50",">50"))
coldat$cc.difference <- coldat$s.score - coldat$g2m.score
rosen_swat@meta.data = coldat

# remove MT and RB high cells to match our own data

dim(rosen_swat)
rosen_swat = subset(rosen_swat, subset = mt.percent < 5.0)
rosen_swat = subset(rosen_swat, subset = ribo.percent < 5.0)
dim(rosen_swat)

DefaultAssay(rosen_swat) = "RNA"
rosen_swat = SplitObject(rosen_swat, split.by = "sample")

# normalise and transform as array job 

njobs=13
interval=ceiling(length(rosen_swat)/njobs)
start=(interval*(tid-1))+1
end=ifelse((interval*tid)<length(rosen_swat), (interval*tid), length(rosen_swat))
rosen_swat=rosen_swat[start:end]

rosen_swat <- lapply(rosen_swat, FUN = ScaleData,  vars.to.regress = c("mt.percent"))
rosen_swat <- lapply(rosen_swat, FUN = SCTransform, vst.flavor = "v2", vars.to.regress = c("mt.percent"))
rosen_swat <- lapply(rosen_swat, FUN = RunPCA)

saveRDS(rosen_swat, paste0(".../rosen_human_sat_split_",tid,".rds"))

"EOF"
