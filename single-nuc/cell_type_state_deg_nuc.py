#!/usr/bin/env python
# coding: utf-8


#import all tools  - if not installed please install them using "pip install toolname --user" or following the instructions on their GitHub page.
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import matplotlib.pyplot as plt
import os
import anndata

import gc
import seaborn as sns



# # Define Required Functions




def filter_by_expression(adata_obj, DE_column, key='rank_genes_groups'):

    def to_tensor(dataframe, columns = [], dtypes = {}, index = False):
        to_records_kwargs = {'index': index}
        if not columns:  # Default to all `dataframe.columns`
            columns = dataframe.columns
        if dtypes:       # Pull in modifications only for dtypes listed in `columns`
            to_records_kwargs['column_dtypes'] = {}
            for column in dtypes.keys():
                if column in columns:
                    to_records_kwargs['column_dtypes'].update({column: dtypes.get(column)})
        return dataframe[columns].to_records(**to_records_kwargs)
    gene_names = pd.DataFrame(adata_obj.uns[key]['names'])
    fraction_in_cluster_matrix = pd.DataFrame(
                np.zeros(gene_names.shape),
                columns=gene_names.columns,
                index=gene_names.index,
            )
    fraction_notin_cluster_matrix = pd.DataFrame(
                np.zeros(gene_names.shape),
                columns=gene_names.columns,
                index=gene_names.index,
            )
    in_cluster_expr_matrix = pd.DataFrame(
                np.zeros(gene_names.shape),
                columns=gene_names.columns,
                index=gene_names.index,
            )
    notin_cluster_matrix = pd.DataFrame(
                np.zeros(gene_names.shape),
                columns=gene_names.columns,
                index=gene_names.index,
            )
    allin_cluster_expr_matrix = pd.DataFrame(
                np.zeros(gene_names.shape),
                columns=gene_names.columns,
                index=gene_names.index,
            )
    allnotin_cluster_matrix = pd.DataFrame(
                np.zeros(gene_names.shape),
                columns=gene_names.columns,
                index=gene_names.index,
            )
    # Create a Adata raw with raw in .X -> We might remove this step - redundant!
    dummy=anndata.AnnData(X=adata_obj.raw.X,
                   var=adata_obj.raw.var,
                   obs=adata_obj.obs)
    # We split the adata by Cell-type. Genes are ordered according to ranks from DEgene analysis
    dummy_split=[dummy[dummy.obs[DE_column].isin([i]),
                  gene_names[i]] for i in gene_names.columns]
    dummy_split_adatasize=[i.shape[0] for i in dummy_split] # How many cells we have in-cluster
    dummy_outcluster=[dummy[~dummy.obs[DE_column].isin([i]),
                       gene_names[i]] for i in gene_names.columns]
    dummy_split_adatasize_outcluster=[i.shape[0] for i in dummy_outcluster] # How many cells we have in out-cluster
    for adata, sizeof_adata, adata_outcluster, sizeof_adata_outcluster in zip(dummy_split, dummy_split_adatasize, 
                                                     dummy_outcluster, dummy_split_adatasize_outcluster):
        cluster_number=[(adata.obs[DE_column].unique()).tolist()][0]

        incluster_counts=np.diff(adata.X.tocsc().indptr)
        percent_expression_incluster=np.array(incluster_counts)/sizeof_adata 
        fraction_in_cluster_matrix[cluster_number[0]]=percent_expression_incluster

        outcluster_counts=np.diff(adata_outcluster.X.tocsc().indptr)
        percent_expression_outcluster=np.array(outcluster_counts)/sizeof_adata_outcluster 
        fraction_notin_cluster_matrix[cluster_number[0]]=percent_expression_outcluster
        # Expression_all in cluster
        allin_cluster_expr_matrix[cluster_number[0]]=np.squeeze(np.asarray(adata.X.mean(axis=0)))
        allnotin_cluster_matrix[cluster_number[0]]=np.squeeze(np.asarray(adata_outcluster.X.mean(axis=0)))
    # Generate the last table results
    result = adata_obj.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    result['pi'] = to_tensor(fraction_in_cluster_matrix)
    result['po'] = to_tensor(fraction_notin_cluster_matrix)
    result['ei'] = to_tensor(allin_cluster_expr_matrix)
    result['eo'] = to_tensor(allnotin_cluster_matrix)
    results =pd.DataFrame(
        {group + '_' + key[:7]: result[key][group]
        for group in groups for key in ['names',"pvals", 'pvals_adj', 'logfoldchanges', 'pi', 'po', "ei", "eo"]})
    return results





path="path_in/"


# # Cell Type DGE




INPUT = path+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad"

adata = sc.read(INPUT)

adata.X=adata.layers["log1p_counts"]

del adata.layers
gc.collect()





adata=adata[~adata.obs.cell_type_am_fine.isin(["Unassigned"])]





adata.obs["cell_type_am"]=adata.obs["cell_type_am"].replace(cell_type_rename)





adata.uns["log1p"]={"base":None}
adata.raw=adata
expression_final=None
for i in adata.obs.cell_type_am_fine.unique():
    tmp=adata.copy()
    tmp.obs["test"]=np.where(tmp.obs.cell_type_am_fine.isin([i]),i,"Everything")
    sc.tl.rank_genes_groups(tmp, "test", method="wilcoxon",use_raw=True,n_genes=36000,reference="Everything")
    expression=filter_by_expression(tmp, 'test', key='rank_genes_groups')
    expression.columns=["Gene","pval","FDR","logFC","PercentageIn","PercentageOut","ExpressionIn","ExpressionOut"]
    expression=expression[["Gene","pval","FDR","logFC","PercentageIn","PercentageOut","ExpressionIn","ExpressionOut"]]
    expression["CellType"]=i
    expression=expression[expression.logFC>0]
    expression=expression[expression.pval<0.05]
    expression=expression.reset_index(drop=True)
    if expression_final is None:
        expression_final=expression.copy()
    else:
        expression_final=pd.concat([expression_final,expression])
    expression.to_csv(path+"/supp_tables/rev/"+i+"_global_markers.csv")

del tmp
expression_final=expression_final.reset_index(drop=True)
expression_final.to_csv(path+"/supp_tables/rev/celltype_global_markers_all.csv")
gc.collect()


# # Cell State DGE




adata.uns["log1p"]={"base":None}
adata.raw=adata
expression_final_all=None
for i in adata.obs.cell_type_am.unique():
    expression_final=None
    if i not in ["Unassigned","Lymphatic","Mast"]: # Lymphatic and Mast do not have cell states. If doing stress vs basal restrict to ["Adipocytes", "APC-ASC", "Endothelial", "Mural"]
        tmp=adata[(adata.obs.cell_type_am.isin([i]))&(~adata.obs.cell_type_am_fine.isin(["Unassigned"]))]
        for j in tmp.obs.cell_state_am.unique():
            tmp.obs["test"]=np.where(tmp.obs.cell_state_am.isin([j]),j,"Everything") # comment to doversus basal
            #tmp.obs["test"]=np.where(tmp.obs.cell_state_am.isin(["AD1","EC1","APC2","Mu1","Mu2"]),"Basal",tmp.obs["cell_state_am"]) # uncomment to do versus basal
            sc.tl.rank_genes_groups(tmp, "test", method="wilcoxon",use_raw=True,n_genes=36000,reference="Everything") # change reference to "Basal" if doing DGE vs basal state
            expression=filter_by_expression(tmp, 'test', key='rank_genes_groups')
            expression.columns=["Gene","pval","FDR","logFC","PercentageIn","PercentageOut","ExpressionIn","ExpressionOut"]
            expression=expression[["Gene","pval","FDR","logFC","PercentageIn","PercentageOut","ExpressionIn","ExpressionOut"]]
            expression["CellType"]=i
            expression["CellState"]=j
            expression=expression[expression.logFC>0]
            expression=expression[expression.pval<0.05]
            if expression_final is None:
                expression_final=expression.copy()
            else:
                expression_final=pd.concat([expression_final,expression])
        expression_final=expression_final.reset_index(drop=True)
        expression_final.to_csv(path+"/supp_tables/rev/"+i+"_state_markers.csv")
    if expression_final_all is None:
        expression_final_all=expression_final.copy()
    else:
        expression_final_all=pd.concat([expression_final_all,expression_final])
del tmp
expression_final_all=expression_final_all.reset_index(drop=True)
expression_final_all.to_csv(path+"/supp_tables/rev/cellstate_markers_all.csv")
gc.collect()







