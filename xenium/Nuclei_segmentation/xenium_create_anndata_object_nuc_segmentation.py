import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import image
from matplotlib import transforms
import seaborn as sns
import numpy as np
import scipy
import gc

import os
import anndata
import logging

import math
import tifffile
import cv2
import sys


path_main="path_in/output-XETG00207__"

# get name from PBS submission request
name = "none"
if("--name" in  sys.argv):
    name = sys.argv[sys.argv.index("--name") + 1]        

   
sample=path_main+name+"/"

# Anndata creation pipeline

transcripts=pd.read_csv(sample+"transcripts.csv.gz") # read xenium transcript table

transcripts=transcripts[transcripts["qv"]>=35] #adjust to desired qv

transcripts=transcripts[~transcripts["feature_name"].str.contains("BLANK|NegControl")] # Remove Blank and negative control probes

transcript_clean=transcripts[(~transcripts["cell_id"].isin(["UNASSIGNED"])) # remove transcripts not assigned to any cells
                             &(transcripts.nucleus_distance<=2)# change distance to be more or less stringent as desired
                            ] 

matrix=pd.crosstab(transcript_clean["cell_id"],transcript_clean["feature_name"])

adata=sc.AnnData(X=np.array(matrix),var=list(matrix.columns),obs=list(matrix.index))

adata.obs["cell_ids"]=adata.obs[0]

del adata.obs[0]

adata.obs.index=adata.obs["cell_ids"]

adata.var["Genes"]=adata.var[0]
del adata.var[0]
adata.var.index=adata.var["Genes"]

adata.X=scipy.sparse.csr_matrix(adata.X)

x_coor=transcript_clean["x_location"].groupby(by=transcript_clean["cell_id"]).mean() # the mean coordinates of all transcripts will be the same as the nuclei centroid due to the short distance.
y_coor=transcript_clean["y_location"].groupby(by=transcript_clean["cell_id"]).mean() # the mean coordinates of all transcripts will be the same as the nuclei centroid due to the short distance.

adata.obs=pd.concat([adata.obs,x_coor],axis=1)
adata.obs=pd.concat([adata.obs,y_coor],axis=1)

obsm=np.stack(adata.obs.apply(lambda row: np.append(int(float(row['x_location'])),int(float(row['y_location']))), axis=1)) # create the array necessary to plot spatial data


gc.collect()


adata.obsm["spatial"]=obsm



#load image files individually
#the image data has been converted to RGB on imageJ. Otherwise it loads as grayscale and plots on a color map
with tifffile.TiffFile(sample+'post_scan_dapi_rgb.tif') as tif:
    post_dapi = tif.asarray()

with tifffile.TiffFile(sample+'post_scan_wga_rgb.tif') as tif:
    post_wga = tif.asarray()
    
with tifffile.TiffFile(sample+'post_scan_dapi_wga_rgb.tif') as tif:
    post_dapi_wga = tif.asarray()


adata.uns["spatial"]=({name: {'images': {'DAPI':post_dapi,
                                        "WGA":post_wga,"DAPI_WGA":post_dapi_wga},
                        'scalefactors':{'tissue_DAPI_scalef': 1.177, # needed to adjust scale as images were re-acquired post-xenium. Default will be 1, but can be verified by plotting afterwards.
                                        'tissue_WGA_scalef': 1.177,
                                        'tissue_DAPI_WGA_scalef': 1.177,
                                        "fiducial_diameter_fullres": 5,'spot_diameter_fullres': 5}}})


adata.obs['n_counts'] = adata.X.sum(axis=1)
adata.obs["n_counts_log"]=np.log1p(adata.obs['n_counts'])

#matrix.replace(0, np.nan, inplace=True) # replace 0 with nan for better plotting down the line
adata.X=np.array(matrix)
adata.raw=adata
sc.pp.log1p(adata)

adata.write(sample+name+"cell_ids_nuc_only_v2.h5ad")



