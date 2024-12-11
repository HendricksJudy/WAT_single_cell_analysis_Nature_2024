#!/usr/bin/env python
# coding: utf-8



import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy



import os
import anndata
import logging

import math

import gc


import warnings
warnings.filterwarnings("ignore")




path_spatial="path_in/"





# load bin object
spatial=sc.read_h5ad(path_spatial+"xenium_bin50_merged_v1.h5ad")





# define stress markers present in spatial dataset
markers=['NAMPT',
 'ETV6',
 'BTG2',
 'MCL1',
 'CCN1',
 'ELOVL5',
 'STAT3',
 'ITGAV',
 'JUNB',
 'ITGB1',
 'JUN',
 'ADAMTS1',
 'OSMR',
 'ATF3',
 'FOSB',
 'ZFP36',
 'HIF1A',
 'RFX2',
 'MYC',
 'CDKN1A',
 'EGFR',
 'EGR1']





# use raw data for scoring
spatial.X=spatial.layers["raw"].copy()
spatial.obs.index=spatial.obs["X_Y_bins"].astype(str)+"__"+spatial.obs["SampleID_simple"].astype(str)#make sure index is unique
spatial.raw=spatial.copy()





# run scoring per condition and merge
final_df=None
for i in spatial.obs.condition.unique():
    tmp_adata=spatial[spatial.obs.condition.isin([i])]
    sc.tl.score_genes(tmp_adata, markers, ctrl_size=50, score_name='Stress_Score_raw', use_raw=None)
    tmp_adata.obs["Stress_Score_quant_con"]=pd.DataFrame(pd.qcut(tmp_adata.obs["Stress_Score_raw"],4,labels=["Q1","Q2","Q3","Q4"])) # get quantiles within condition
    if final_df is None:
        final_df=tmp_adata.obs[["Stress_Score_raw","Stress_Score_quant_con"]]
    else:
        final_df=pd.concat([final_df,tmp_adata.obs[["Stress_Score_raw","Stress_Score_quant_con"]]],axis=0)





# reindex merged datframe to match spatial.obs
final_df=final_df.reindex(spatial.obs.index)





spatial.obs["Stress_Score_raw"]=final_df["Stress_Score_raw"]





spatial.obs["Stress_Score_quant_con"]=final_df["Stress_Score_quant_con"]





# calculate quantiles across all conditions
spatial.obs["Stress_Score_quant_all"]=pd.DataFrame(pd.qcut(spatial.obs["Stress_Score_raw"],4,labels=["Q1","Q2","Q3","Q4"]))













