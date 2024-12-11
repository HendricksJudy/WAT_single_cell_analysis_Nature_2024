#!/usr/bin/env python
# coding: utf-8


import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import image
from matplotlib import transforms
import seaborn as sns
import numpy as np
import scipy
import bbknn as bk

import os
import anndata
import logging

import math
import tifffile
import cv2
import gc
import harmonypy as hm





sc.settings.set_figure_params(dpi=300,fontsize=10)





path_main="/rds/general/user/aalmeid2/projects/lms-scott-raw/live/Xenium/wl_project/"





# load global object
adata=sc.read_h5ad(path_main+"xenium_merged_global_v3.h5ad")





# get all files with cell state annotation
files=[x for x in os.listdir() if "xenium_" in x and "cell_state_v3.csv" in x]





# load and merge all files
final_df=None
for i in files:
    tmp=pd.read_csv(i,index_col=0)
    if final_df is None:
        final_df=tmp.copy()
    else:
        final_df=pd.concat([final_df,tmp])
    del tmp
    gc.collect()





#concatenate cell type annotation from adata.obs
#this also sorts the indexes into the correct order
df=pd.concat([adata.obs["cell_type"],final_df],axis=1)





# label Mast cells in the cell state column
df["cell_state"]=np.where(df["cell_type"].isin(["Mast"]),"Mast",df["cell_state"])





# add cell state column to adata.obs
adata.obs=pd.concat([adata.obs,df["cell_state_simple"]],axis=1)





# add Condition column to adata.obs
sample_dict={
"A1":"Obese",
"A5":"Weightloss",
"A11":"Lean",
"A4":"Obese",
"A8":"Weightloss",
"A16":"Lean",
"A2":"Obese",
"A3":"Obese",
"A6":"Weightloss",
"A7":"Weightloss",
"A13":"Lean",
"A14":"Lean"
}
adata.obs["condition"]=adata.obs["SampleID_simple"].replace(sample_dict)





#save object 
adata.write(path_main+"xenium_merged_global_v3.h5ad")







