import scanpy as sc
import pandas as pd
import numpy as np
import scipy
from scipy import spatial
import os
import anndata



import sys

path_main="path_in/"

adata=sc.read_h5ad(path_main+"xenium_merged_global_v3.h5ad") # load anndata object will all samples

adata.obs["x_y"]=adata.obs["x_location"].astype(str)+"__"+adata.obs["y_location"].astype(str)

for j in adata.obs["SampleID_simple"].unique(): # query needs to be per sample
    final_distance=[]
    final_order=[]
    final_query=[]
    final_sample=[]
    adata_tmp=adata[adata.obs.SampleID_simple.isin([j])]
    test=np.column_stack((adata_tmp.obs.x_location,adata_tmp.obs.y_location))
    tree = spatial.KDTree(test)
    for i in adata_tmp.obs["x_y"]:
        x,y=i.split("__")
        x=float(x)
        y=float(y)
        qtest=tree.query([(x,y)],k=len(adata_tmp))
        distance=np.squeeze(qtest[0]) #get all the KD distances
        order=np.squeeze(qtest[1]) #get all the target cells
        query=np.full(shape=len(distance),fill_value=i)
        sample=np.full(shape=len(distance),fill_value=j)
        final_distance=np.append(final_distance,distance)
        final_order=np.append(final_order,order)
        final_query=np.append(final_query,query)
        final_sample=np.append(final_sample,sample)
    df=pd.DataFrame([final_distance,final_order,final_query,final_sample])
    df=df.T
    df.columns=["distance","target","query","sample"]
    cells=pd.Series(adata_tmp.obs["cell_ids"].values,index=adata_tmp.obs["x_y"])
    df["query"]=df["query"].replace(cells)
    order=pd.Series(adata_tmp.obs["cell_ids"].values).to_dict()
    df["target"]=df["target"].replace(order)
    df.to_csv(path_main+j+"_all_distance.csv") # export master table with distances between all cells within samples. Cell identity (cell type, cell state, etc.) can be changed afterwards depending on needs