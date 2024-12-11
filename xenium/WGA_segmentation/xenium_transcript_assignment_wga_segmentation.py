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
from scipy import spatial

import os
import anndata
import logging

import math
import tifffile
import cv2
import gc
import sys



samples=[
"0014085__A1__20240208__122452",
"0014085__A5__20240208__122452",
"0014085__A11__20240208__122452",
"0014253__A4__20240208__122452",
"0014253__A8__20240208__122452",
"0014253__A16__20240208__122452",
"0014271__A2__20240222__114657",
"0014271__A6__20240222__114657",
"0014271__A13__20240222__114657",
"0014373__A3__20240222__114657",
"0014373__A7__20240222__114657",
"0014373__A14__20240222__114657",
]




path_main="path_in/"




for name in samples:
    id1,id2,id3,id4=name.split("__")

    folder="/output-XETG00207__"

    sample=path_main+folder+name+"/"
    transcripts=pd.read_csv(sample+"transcripts.csv.gz")# load xenium transcript table

    transcripts=transcripts[transcripts["qv"]>=35] #set desired quality treshold

    transcripts=transcripts[~transcripts["feature_name"].str.contains("BLANK|NegControl")] # remove Blank and Negative Control probes

    transcripts["x_location"]=transcripts["x_location"]*1.177 #adjut coordinates due to post-xenium scan
    transcripts["y_location"]=transcripts["y_location"]*1.177 #adjut coordinates due to post-xenium scan

    adata=sc.read_h5ad(path_main+"xenium_merged_global_loose.h5ad") # load a nuclei segmentation object with cell type clustering. This will be used to remove any non-adipocyte transcripts
    
    adata=adata[(adata.obs.SampleID_simple.isin([id2]))&(~adata.obs.cell_type.isin(["Adipocytes"]))] # subset to non-adipocytes and to the sample of interest

    adip_seg=pd.read_csv(path_main+"wga_segmentation/"+id2+"_ij_result.csv") #get table with WGA segmentation pixel coordinates from imageJ
    seg_metrics=pd.read_csv(path_main+"wga_segmentation/"+id2+"_ij_result_metrics.csv") # get table with "measure" results from imageJ

    seg_metrics=seg_metrics[(seg_metrics.Area>=1000)&(seg_metrics.Area<=25000)] # remove objects that are too small or too big
    adip_seg=adip_seg[adip_seg.object.isin(seg_metrics.iloc[:,0])]
    adip_seg=adip_seg.reset_index(drop=True)
    non_adip=adata.obs.cell_ids.unique()
    del adata
    gc.collect()
    transcripts_non_adip=transcripts[(transcripts.cell_id.isin(non_adip))&(transcripts.nucleus_distance<=10)] # adjust distance to remove more or fewer genes. The smaller the distance, the higher the chance of including noise
    transcripts_adip=transcripts[~transcripts.transcript_id.isin(transcripts_non_adip.transcript_id.unique())]
    transcripts_non_adip=transcripts_non_adip[(transcripts_non_adip.nucleus_distance<=2)] # subset the non-adip transcripts to the final desired distance to create the mixed object at the end
    
    # assign transcripts to segmented adipocytes
    
    test=np.column_stack((adip_seg.xpoints,adip_seg.ypoints))
    tree = spatial.KDTree(test)

    del transcripts
    gc.collect()
    transcripts_adip["x_y"]=transcripts_adip["x_location"].astype(str)+"__"+transcripts_adip["y_location"].astype(str)
    test=transcripts_adip.apply(lambda row: np.squeeze(tree.query([(row["x_location"],row["y_location"])],k=10)[0]), axis=1) #k=10 assignes to the 10 nearest objects. Used to save computational time. Increase to get more.
    test1=transcripts_adip.apply(lambda row: np.squeeze(tree.query([(row["x_location"],row["y_location"])],k=10)[1]), axis=1) #k=10 assignes to the 10 nearest objects. Used to save computational time. Increase to get more.
    test2=transcripts_adip.apply(lambda row: np.full(shape=len(np.squeeze(tree.query([(row["x_location"],row["y_location"])],k=10)[1])) #k=10 assignes to the 10 nearest objects. Used to save computational time. Increase to get more.
                                                     ,fill_value=row["x_y"]), axis=1)
    test3=transcripts_adip.apply(lambda row: np.full(shape=len(np.squeeze(tree.query([(row["x_location"],row["y_location"])],k=10)[1])) #k=10 assignes to the 10 nearest objects. Used to save computational time. Increase to get more.
                                                     ,fill_value=row["feature_name"]), axis=1)
    distance_array=np.concatenate(np.asarray(test))
    order_array=np.concatenate(np.asarray(test1))
    transcript_array=np.concatenate(np.asarray(test2))
    feature_array=np.concatenate(np.asarray(test3))
    tmp_df=pd.DataFrame([distance_array,order_array,transcript_array,feature_array]).T
    tmp_df.columns=["nucleus_distance","cell_id","x_y","feature_name"] # called nucleus_distance just to merge with other transcripts afterwards. However it represents distance to nearest pixel.
    obj_dict=pd.Series(adip_seg["object"].values,index=adip_seg.index).to_dict() # create a dictionary with object name
    tmp_df["cell_id"]=[obj_dict[x] for x in tmp_df["cell_id"]]
    tmp_df["cxy"]=tmp_df["cell_id"].astype(str)+"__"+tmp_df["x_y"].astype(str)
    tmp_df=tmp_df[tmp_df.nucleus_distance<=2] # subset to the desired maximum distance between trancript and segmented pixel
    tmp_df=tmp_df.drop_duplicates("cxy") # remove any transcripts that have been assigned more than once to the same object
    trans_dict=pd.Series(transcripts_adip.feature_name.values,index=transcripts_adip.x_y).to_dict()
    tmp_df[["x_location","y_location"]]=tmp_df.x_y.str.split("__",expand=True)

    transcripts=pd.concat([transcripts_non_adip,tmp_df],axis=0) # merger non-adipocyte and adipocyte transcript tables
    
    area_dict=pd.Series(seg_metrics["Area"].values,index=seg_metrics.iloc[:,0].values).to_dict() # add area
    x_dict=pd.Series(seg_metrics["X"].values,index=seg_metrics.iloc[:,0].values).to_dict() # add segmented adipocyte centroid
    y_dict=pd.Series(seg_metrics["Y"].values,index=seg_metrics.iloc[:,0].values).to_dict() # add segmented adipocyte centroid
    # add other desired morphological measurements
    circ_dict=pd.Series(seg_metrics["Circ."].values,index=seg_metrics.iloc[:,0].values).to_dict()
    round_dict=pd.Series(seg_metrics["Round"].values,index=seg_metrics.iloc[:,0].values).to_dict()
    feret_dict=pd.Series(seg_metrics["Feret"].values,index=seg_metrics.iloc[:,0].values).to_dict()
    
    # create anndata object
    matrix=pd.crosstab(transcripts["cell_id"],transcripts["feature_name"])

    adata=sc.AnnData(X=np.array(matrix),var=list(matrix.columns),obs=list(matrix.index))

    adata.obs["cell_ids"]=adata.obs[0]

    del adata.obs[0]

    adata.obs.index=adata.obs["cell_ids"]

    adata.var["Genes"]=adata.var[0]
    del adata.var[0]
    adata.var.index=adata.var["Genes"]

    adata.X=scipy.sparse.csr_matrix(adata.X)
    transcripts["x_location_adip"]=transcripts["cell_id"].replace(x_dict)
    transcripts["y_location_adip"]=transcripts["cell_id"].replace(y_dict)

    transcripts["x_location_final"]=np.where(transcripts.cell_id.isin(transcripts_non_adip.cell_id.unique()),transcripts["x_location"],transcripts["x_location_adip"])
    transcripts["y_location_final"]=np.where(transcripts.cell_id.isin(transcripts_non_adip.cell_id.unique()),transcripts["y_location"],transcripts["y_location_adip"])

    transcripts["x_location_final"]=transcripts["x_location_final"].astype(float)
    transcripts["y_location_final"]=transcripts["y_location_final"].astype(float)

    x_coor=transcripts["x_location_final"].groupby(by=transcripts["cell_id"]).mean() # the mean coordinates of all transcripts will be the same as the nuclei centroid due to the short distance.
    y_coor=transcripts["y_location_final"].groupby(by=transcripts["cell_id"]).mean() # the mean coordinates of all transcripts will be the same as the nuclei centroid due to the short distance.

    adata.obs=pd.concat([adata.obs,x_coor],axis=1)
    adata.obs=pd.concat([adata.obs,y_coor],axis=1)

    adata.obs.index=id2+"_"+adata.obs["cell_ids"].astype(str)
    obsm=np.stack(adata.obs.apply(lambda row: np.append(int(row['x_location_final']),int(row['y_location_final'])), axis=1))

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
                                            "fiducial_diameter_fullres": 5,'spot_diameter_fullres': 5}}})
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.obs["n_counts_log"]=np.log1p(adata.obs['n_counts'])

    adata.obs["Area"]=adata.obs["cell_ids"].replace(area_dict)
    adata.obs["Circularity"]=adata.obs["cell_ids"].replace(circ_dict)
    adata.obs["Roundness"]=adata.obs["cell_ids"].replace(round_dict)
    adata.obs["Feret"]=adata.obs["cell_ids"].replace(feret_dict)
    adata.obs["cell_ids"]=adata.obs.index

    adata.X=np.array(matrix)
    adata.raw=adata
    sc.pp.log1p(adata)
    # transforming measurements into strings because non-adipocytes are NaN which can affect how anndata is saved
    adata.obs["Area"]="s_"+adata.obs["Area"].astype(str) 
    adata.obs["Circularity"]="s_"+adata.obs["Circularity"].astype(str)
    adata.obs["Roundness"]="s_"+adata.obs["Roundness"].astype(str) 
    adata.obs["Feret"]="s_"+adata.obs["Feret"].astype(str)
    adata.write(sample+"adip_seg.h5ad")
    gc.collect()
gc.collect()







