import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import image
from matplotlib import transforms
import seaborn as sns
import numpy as np
import scipy
from scipy import spatial

import os
import anndata
import logging

import math
import tifffile
import cv2
import gc
import sys


path_main="path_in/"

# get sample name from PBS submission
name = "none"
if("--name" in  sys.argv):
    name = sys.argv[sys.argv.index("--name") + 1]        
    
id1,id2,id3,id4=name.split("__")


folder="/output-XETG00207__"
    


sample=path_main+folder+name+"/"

transcripts=pd.read_csv(sample+"transcripts.csv.gz") # load xenium transcript table

transcripts=transcripts[transcripts["qv"]>=35] #set desired quality treshold

transcripts=transcripts[~transcripts["feature_name"].str.contains("BLANK|NegControl")] # remove Blank and Negative Control probes

transcripts["x_location"]=transcripts["x_location"]*1.177 #adjut coordinates due to post-xenium scan
transcripts["y_location"]=transcripts["y_location"]*1.177 #adjut coordinates due to post-xenium scan


raw=pd.DataFrame()
raw["X"]=transcripts["x_location"]
raw["Y"]=transcripts["y_location"]
raw["Z"]=transcripts["z_location"]
raw["Gene"]=transcripts["feature_name"]

n_bins=50 # choose desired Bin dimension (ensure pixel dimensions are 1 micron, otherwise adjust value accordingly)
raw["X_bin"]=(round(raw["X"] / n_bins)) * n_bins
raw["Y_bin"]=(round(raw["Y"] / n_bins)) * n_bins
    
raw["X_Y_bins"]=raw["X_bin"].astype(str)+"_"+raw["Y_bin"].astype(str)

raw["nucleus"]=transcripts["overlaps_nucleus"]


matrix=pd.crosstab(raw["X_Y_bins"],raw["Gene"])

adata=sc.AnnData(X=np.array(matrix),var=list(matrix.columns),obs=list(matrix.index))

adata.obs["X_Y_bins"]=adata.obs[0]

del adata.obs[0]

adata.obs.index=adata.obs["X_Y_bins"]
adata.obs[['X_coor', 'Y_coor']] = adata.obs["X_Y_bins"].apply(lambda x: pd.Series(str(x).split("_")))

adata.var["Genes"]=adata.var[0]
del adata.var[0]
adata.var.index=adata.var["Genes"]

adata.X=scipy.sparse.csr_matrix(adata.X)

obsm=np.stack(adata.obs.apply(lambda row: np.append(int(float(row['X_coor'])),int(float(row['Y_coor']))), axis=1)) # create array for plotting based on bin centroid

gc.collect()

adata.obsm["spatial"]=obsm
# load post-xenium images
#the image data has been converted to RGB on imageJ. Otherwise it loads as grayscale and plots on a color map

with tifffile.TiffFile(sample+'post_scan_dapi_rgb.tif') as tif:
    post_dapi = tif.asarray()

with tifffile.TiffFile(sample+'post_scan_wga_rgb.tif') as tif:
    post_wga = tif.asarray()
    
with tifffile.TiffFile(sample+'post_scan_dapi_wga_rgb.tif') as tif:
    post_dapi_wga = tif.asarray()
    
with tifffile.TiffFile(path_main+"wga_segmentation/"+id2+'_binary_rgb.tif') as tif:
    post_wga_seg = tif.asarray()


adata.uns["spatial"]=({name: {'images': {'DAPI':post_dapi,
                                        "WGA":post_wga,"DAPI_WGA":post_dapi_wga,"WGA_seg":post_wga_seg},
                        'scalefactors':{'tissue_DAPI_scalef': 1,'tissue_WGA_seg_scalef': 1,
                                        'tissue_WGA_scalef': 1,'tissue_DAPI_WGA_scalef': 1,
                                        "fiducial_diameter_fullres": 50,'spot_diameter_fullres': 50}}})

adata.obs['n_counts'] = adata.X.sum(axis=1)
adata.obs["n_counts_log"]=np.log1p(adata.obs['n_counts'])


#matrix.replace(0, np.nan, inplace=True) # replace 0 with nan for better plotting down the line (optional - some applications it needs to be reverted if used)
adata.X=np.array(matrix)
adata.raw=adata



adata.write(sample+"adip_bin50_v1.h5ad")
