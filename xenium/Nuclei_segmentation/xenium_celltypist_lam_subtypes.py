#!/usr/bin/env python
# coding: utf-8



import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import celltypist as ct


import gc





sc.settings.set_figure_params(dpi=300,fontsize=10,dpi_save=300)





path="/path_in/"





# load spatial data
spatial=sc.read_h5ad(path_spatial+"xenium_merged_global_v3.h5ad")
spatial.X=spatial.layers["raw"]





# normalize and log transform spatial data - this is required by celltypist
sc.pp.normalize_total(spatial, target_sum=1e4)
sc.pp.log1p(spatial)





spatial.raw=spatial.copy()





# subset only to LAM cluster
# "best match" will be used for the prediction. This step prevents other cell-types/states from being used for the prediction.
lams=spatial[spatial.obs.cell_state_simple.isin(["MYE2"])]





lams.raw=lams.copy()





# run prediction using the 2nd round model
predictions = ct.annotate(lams, model = "human_adip_lams_celltypist_fine_SGD_alpha001_r2_prob05_01_spatial.pkl", majority_voting = True, mode = 'best match')
spatial = predictions.to_adata(insert_prob = False)
gc.collect()







