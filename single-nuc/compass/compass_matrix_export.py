#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import scanpy as sc
import gc

from collections import OrderedDict

import re
import scipy
import scipy.sparse

import anndata


# In[2]:


def get_df(scanpy_object, column_tosplitby, sum_or_mean):
    d = {}
    for cluster_number in np.unique(scanpy_object.obs[column_tosplitby].values):
        scanpy_object_subset = scanpy_object[scanpy_object.obs[column_tosplitby].isin([cluster_number])]

        if sum_or_mean=="mean":
            d[cluster_number] = np.squeeze(np.asarray(scanpy_object_subset.raw.X.mean(axis=0))) 
        elif sum_or_mean=="sum":
            d[cluster_number] = np.squeeze(np.asarray(scanpy_object_subset.raw.X.sum(axis=0)))
        elif sum_or_mean=="pct":
            d[cluster_number] = np.squeeze(np.asarray((scanpy_object_subset.raw.X > 0).mean(axis=0)))   
        del scanpy_object_subset
        
    return d


# In[3]:


path_out="/path_out/compass/"


# In[4]:


path_in="/path_in/"
INPUT = path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad"


# In[5]:


adata = sc.read(INPUT)


# In[6]:


adata


# In[7]:


adata.X.max()


# In[8]:


adata.obs["pool_key"]=adata.obs["sample"] # get expression per sample


# In[9]:


adata.obs["pool_key"].value_counts()


# In[10]:


adata.X=adata.layers["log1p_counts"]


# In[11]:


adata.X.max()


# In[12]:


adata.X


# In[13]:


adata_orig=adata[adata.obs.cell_type_am_fine.isin(["Adipocytes"])] # change cell_type accordingly
adata_orig.raw=adata_orig
x=get_df(adata_orig, 'pool_key', 'mean') # get mean expression matrix per sample
x=pd.DataFrame(x)
x.index=adata_orig.var.index


# In[14]:


x


# In[15]:


x.to_csv(path_out+"log1p_mean_adipocytes/expression_matrix_log1p_per_sample_mean.tsv",sep="\t")


# In[ ]:




