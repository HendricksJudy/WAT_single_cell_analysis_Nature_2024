#!/usr/bin/env python
# coding: utf-8



import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
import anndata2ri
import rpy2.rinterface_lib.callbacks
import logging
import bbknn as bk
import gc





# Ignore R warning messages
#Note: this can be commented out to get more verbose R output
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)





anndata2ri.activate()





get_ipython().run_line_magic('load_ext', 'rpy2.ipython')





get_ipython().run_cell_magic('R', '', 'library(miloR)\nlibrary(SingleCellExperiment)\nlibrary(scater)\nlibrary(scran)\nlibrary(dplyr)\nlibrary(patchwork)\nlibrary(igraph)')





path="path_in/"
#INPUT = path+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad" # use for cell_type_fine analysis
INPUT = path+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_Adipocytes_noDoublets.h5ad" # change to the desired cell type





adata = sc.read(INPUT)



adata.X=adata.layers["log1p_counts"]
del adata.layers
adata.raw=adata
gc.collect()


# # Add missing co-variates




covariate=pd.read_csv("phe_snRNAseq_WAT_phase_4_incomplete.csv")





df=pd.DataFrame([adata.obs["sample"],adata.obs["age"],adata.obs["sex"],
                 adata.obs["race"],adata.obs["ethnicity"],adata.obs["bmi"]]).T





df["barcode"]=df.index





df.index=df["sample"]





covariate.index=covariate["SID"]





df["age"][df["sample"].isin(covariate["SID"])]=df["sample"]
df["ethnicity"][df["sample"].isin(covariate["SID"])]=df["sample"]
df["sex"][df["sample"].isin(covariate["SID"])]=df["sample"]
df["bmi"][df["sample"].isin(covariate["SID"])]=df["sample"]





df





age_rename=zip(covariate.SID,covariate.Age)
age_rename=dict(age_rename)
sex_rename=zip(covariate.SID,covariate.Sex)
sex_rename=dict(sex_rename)
eth_rename=zip(covariate.SID,covariate.Ethnicity)
eth_rename=dict(eth_rename)
bmi_rename=zip(covariate.SID,covariate.BMI)
bmi_rename=dict(bmi_rename)





df["age"]=df["age"].replace(age_rename)
df["sex"]=df["sex"].replace(sex_rename)
df["ethnicity"]=df["ethnicity"].replace(eth_rename)
df["bmi"]=df["bmi"].replace(bmi_rename)





df["ethnicity"]=df["ethnicity"].replace({
    "Caucasian":"EW",
    "Black":"AF"
})





adata.obs["ethnicity2"]=np.asarray(df["ethnicity"]).astype(str)
adata.obs["age2"]=np.asarray(df["age"]).astype(float)
adata.obs["sex2"]=np.asarray(df["sex"]).astype(str)
adata.obs["bmi2"]=np.asarray(df["bmi"]).astype(float)





df["FG"]=np.nan
df["FG"][df["sample"].isin(covariate["SID"])]=df["sample"]
fg_rename=zip(covariate.SID,covariate.FG)
fg_rename=dict(fg_rename)
df["FG"]=df["FG"].replace(fg_rename)
adata.obs["FG"]=np.asarray(df["FG"]).astype(float)

df["FI"]=np.nan
df["FI"][df["sample"].isin(covariate["SID"])]=df["sample"]
fg_rename=zip(covariate.SID,covariate.FI)
fg_rename=dict(fg_rename)
df["FI"]=df["FI"].replace(fg_rename)
adata.obs["FI"]=np.asarray(df["FI"]).astype(float)


# # Milo Pipeline - Obese vs Lean




adata_sub=adata[~adata.obs.condition2.isin(["Weightloss"])] # for Weightloss remove Lean instead. To analyse within a single condition, subset specifically to the desired condition.
gc.collect()





adata_sub=adata_sub[(adata_sub.obs.dataset.isin(["scott"]))&(~adata_sub.obs.cell_state_am.isin(["Unassigned"]))]





bk.bbknn(adata_sub,n_pcs=40,local_connectivity=1,batch_key="sample",neighbors_within_batch=1)





adata_no_knn = adata_sub.copy()
adata_no_knn.obsp = None
adata_no_knn.uns.pop("neighbors")
adata_no_knn





get_ipython().run_cell_magic('R', '-i adata_no_knn', 'adata_no_knn')





get_ipython().run_cell_magic('R', '', 'milo <- Milo(adata_no_knn)\nmilo')





## Save the binary connectivity matrix
knn_adjacency = adata_sub.obsp["connectivities"]





knn_adjacency






get_ipython().run_cell_magic('R', '-i knn_adjacency', '\nmilo_graph <- buildFromAdjacency(knn_adjacency, k=20, is.binary=TRUE)\ngraph(milo) <- miloR::graph(milo_graph)')





get_ipython().run_cell_magic('R', '', 'milo <- buildGraph(milo, k=20, d=30)')





design_df = adata_sub.obs[["sample", "condition2"]] # add any other desired co-variates
design_df.drop_duplicates(inplace=True)
design_df.index = design_df['sample']
design_df





get_ipython().run_cell_magic('R', '-i design_df -o DA_results', '## Define neighbourhoods\nmilo <- makeNhoods(milo, prop = 0.1, k = 20, d=30, refined = TRUE)\n\n## Count cells in neighbourhoods\nmilo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="sample")\n\n## Calculate distances between cells in neighbourhoods\n## for spatial FDR correction\nmilo <- calcNhoodDistance(milo, d=30)\n## Test for differential abundance\nDA_results <- testNhoods(milo, \n                         design = ~ condition2, #change to desired co-variates. When doing multiple, the p-value provided is for the right-most variable. For example, bmi2+FI represents FI adjusted for BMI\n                         design.df = design_df)')





gc.collect()





plt.plot(DA_results.logFC, -np.log10(DA_results.SpatialFDR), '.');
plt.xlabel("log-Fold Change");
plt.ylabel("- log10(Spatial FDR)")





get_ipython().run_cell_magic('R', '', 'milo <- buildNhoodGraph(milo)')





get_ipython().run_cell_magic('R', '-w 1000 -h 800', 'plotNhoodGraphDA(milo, DA_results, alpha=0.05)')





get_ipython().run_cell_magic('R', '-oc', 'DA_results <- annotateNhoods(milo, DA_results, coldata_col = "cell_state_am")')





get_ipython().run_cell_magic('R', '', 'ggplot(DA_results, aes(cell_state_am_fraction)) + geom_histogram(bins=50)')





get_ipython().run_cell_magic('R', '-o DA_results', 'DA_results$cell_state_am <- ifelse(DA_results$cell_state_am < 0.9, "Mixed", DA_results$cell_state_am)')





DA_results





DA_results.to_csv("adip_milo_condition_scott_only_lean.csv")





get_ipython().run_cell_magic('R', '', 'plotDAbeeswarm(DA_results, group.by = "cell_state_am")')





del adata_sub
gc.collect()

