#!/usr/bin/env python
# coding: utf-8



import pandas as pd
import numpy as np
import gc





path_in="input_path/"

path_out="output_path/"





# we have 10 adjacency outputs from pyScenic (from 10 subsampled obects per analysis)
# we create a single adjacency file with the mean across all 10 runs
# load and merge all runs
final_adj=pd.DataFrame()
for i in range(1,11):
    tmp=pd.read_csv(path_in+"global_all_stringent_sub_compass_adipocytes"+str(i)+"adj.tsv",sep="\t")
    tmp["TF_target"]=tmp["TF"]+"_"+tmp["target"]
    tmp.index=tmp["TF_target"]
    del tmp["TF_target"]
    del tmp["TF"]
    del tmp["target"]
    tmp=tmp.add_suffix("_sub"+str(i))
    final_adj=pd.concat([final_adj,tmp],axis=1)





counts=final_adj.T.count()





# assign 0 to runs were a pair was not correlated
final_adj=final_adj.fillna(0)





# calculate the mean adjacency score
means=final_adj.T.mean()





final_adj_single=pd.DataFrame(means.index)





final_adj_single





#split TF and Target into different columns
final_adj_single[["TF","target"]]=final_adj_single["TF_target"].str.split('_',expand=True)





final_adj_single["importance"]=np.asarray(means)





final_adj_single["counts"]=np.asarray(counts)





# Keep only pairs that were found in at least 2 runs
final_adj_single=final_adj_single[final_adj_single.counts>=2]





del final_adj_single["TF_target"]
del final_adj_single["counts"]





final_adj_single=final_adj_single.reset_index(drop=True)





final_adj_single.to_csv(path_out+"adj_merged_compass_adipocytes.tsv",sep="\t",index=False)







