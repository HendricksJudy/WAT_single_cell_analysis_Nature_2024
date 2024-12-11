#!/usr/bin/env python
# coding: utf-8



import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import matplotlib.pyplot as plt
import os
import bbknn as bk
import anndata 
import harmonypy as hm
import gc
get_ipython().run_line_magic('matplotlib', 'inline')





sc.settings.set_figure_params(dpi=300,fontsize=10,dpi_save=300)


# ### Define functions to get DGE tables from anndata




def makevector_topmarkers(anndata, top_n, ldFC_cutoff=0.5, pvals_adj_cutoff=10**-30, expression_in=1,expression_out=1,percentage_out=1,percentage_in=0): 
    results=anndata.uns['rank_genes_groups']
    clusters = results['names'].dtype.names
    top_genes=[]
    for cluster in clusters:
        topgenes_cluster=results['names'][cluster][np.where((results['logfoldchanges'][cluster]>=ldFC_cutoff) & 
                                                            (results['pvals_adj'][cluster]<(pvals_adj_cutoff)) &
                                                            (results['eo'][cluster]<=expression_out) &
                                                            (results['po'][cluster]<=percentage_out) &
                                                            (results['pi'][cluster]>=percentage_in) &
                                                           (results['ei'][cluster]>=expression_in))][0:top_n]
        top_genes.extend(topgenes_cluster)
    return(top_genes)





def filter_by_expression(adata_obj, DE_column, key='rank_genes_groups'):
#adata_obj=fbs_new_harmony.copy()
#key='rank_genes_groups'
#DE_column='leiden'
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
        {group + '_' + key[:2]: result[key][group]
        for group in groups for key in ['names', 'pvals_adj', 'logfoldchanges', 'pi', 'po', "ei", "eo"]})
    return results


# # Load and prepare data




path_in="path_in/"





seed=5





adata=sc.read_h5ad(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna.h5ad")





gc.collect()





adata





#verify if adata.X is log transformed
adata.X.max()





#verify if adata.raw.X is log transformed
adata.raw.X.max()





# get raw matrix to verify data normalization
matrix=pd.DataFrame(adata.raw.X.todense(),columns=adata.var.index,index=adata.obs.index) 





matrix.T.sum()





del matrix
gc.collect()





# get raw log to verify data normalization

matrix=pd.DataFrame(adata.X.todense(),columns=adata.var.index,index=adata.obs.index) 





matrix.T.sum()


# ### Normalize data




# store raw counts in a anndata layer
adata.layers["raw_counts"]=adata.raw.X





adata.X=adata.layers["raw_counts"].copy()





# normalize to 10000 counts per cells
sc.pp.normalize_total(adata, target_sum=1e4)





# get normalized matrix to verify normalizations
matrix=pd.DataFrame(adata.X.todense(),columns=adata.var.index,index=adata.obs.index) 





matrix.T.sum()





# add normalzied counts to layers
adata.layers["norm_counts"]=adata.X.copy()





#log transform the data
sc.pp.log1p(adata)





# add normalzied log counts to layers
adata.layers["log1p_counts"]=adata.X.copy()





adata.raw=adata.copy()





gc.collect()


# # Integration and Clustering




#regress effects from mitochondrial, ribossomal and total counts
sc.pp.regress_out(adata,["mt.percent","nCount_RNA","ribo.percent"])





sc.pp.highly_variable_genes(adata)





sc.tl.pca(adata, svd_solver='arpack',n_comps=40,random_state=seed,use_highly_variable=True)


# ### Run Harmony to correct PCA




pca = adata.obsm['X_pca']
batch = adata.obs['sample']
meta_data = adata.obs

ho = hm.run_harmony(pca, meta_data, ['sample'], epsilon_harmony = -float('Inf'), max_iter_kmeans = 15, tau=5, max_iter_harmony = 25)





res = pd.DataFrame(ho.Z_corr)
res=res.T
adata.obsm['X_pca'] = res.values


# ### Run bbknn to get batch corrected neighbors




bk.bbknn(adata,n_pcs=40,local_connectivity=1,batch_key="sample",neighbors_within_batch=2)


# # Broad Cell Type Clustering and visualization




sc.tl.umap(adata,spread=1,min_dist=0.15)





# low resolution to identify only major cell types
sc.tl.leiden(adata, key_added="cell_type_am",resolution=0.15,random_state=seed,n_iterations=-1)





sc.pl.umap(adata,color="cell_type_am")





sc.pl.umap(adata,color="sample")





sc.pl.umap(adata,color=["dataset","condition"])
# Rosen's data does not contain explicit information on condition. Will be added in downstream analysis based on BMI





sc.pl.umap(adata,color=["nCount_RNA","mt.percent","ribo.percent"])


# # Cluster DGE




sc.tl.rank_genes_groups(adata, "cell_type_am", method="t-test_overestim_var",use_raw=True)
sc.tl.dendrogram(adata,groupby="cell_type_am")





sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, groupby="cell_type_am",mean_only_expressed=True,standard_scale="var")





#this adds the expression in and expression out information to adata and creates a table that can be exported

expression=filter_by_expression(adata, 'cell_type_am', key='rank_genes_groups')





expression





#Dotplot with filtering options for unbiased markers
#Arguments for out cluster are <=,i.e, a percent_out=1 means no filtering will occur based on percentage outside the cluster
#Arguments for in cluster are >=,i.e, a percent_in=0 means no filtering will occur based on percentage inside the cluster


sc.pl.dotplot(adata,
              makevector_topmarkers(adata, 
                                    5, ldFC_cutoff=1, pvals_adj_cutoff=10**-10, expression_in=0.5,expression_out=0.3,percentage_out=0.3,percentage_in=0.5),
              groupby='cell_type_am', use_raw=True,standard_scale="var")





# add regressed counts to layers
adata.layers["regressed_log1p_counts"]=adata.X





#save
adata.write(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated.h5ad")


# # Annotate Clusters




#hand picked markers
marker_genes={"Adipocytes":["ADIPOQ","GPAM"],
              "Endothelial":["VWF","CDH5"],
              "Lymphatic":["PROX1"],
              "APC":["DCN","PDGFRA"],
              "Immune":["PTPRC"],
              "Myeloid":["MRC1"],
              "Lymphoid":["CD247"],
              "Mast":["KIT","CPA3"],
              "Smooth-Musle":["MYH11"],
              "Pericytes":["KCNJ8"]
             }





sc.pl.dotplot(adata,
              marker_genes,
              groupby='cell_type_am', use_raw=True,standard_scale="var")





#rename clusters to respective cell type
adata.obs.cell_type_am=adata.obs.cell_type_am.replace(
    {"0":"Adipocytes",
     "1":"APC-ASC",
     "2":"Myeloid",
     "3":"Endothelial",
     "4":"Lymphoid",
     "5":"Mural",
     "6":"Mast",
     "7":"Lymphatic"
    })





marker_genes=["ADIPOQ","GPAM","DCN","PDGFRA","PTPRC","MRC1","CD68","VWF","CDH5","CD247","NKG7","MYH11","KCNJ8","KIT","CPA3","PROX1","CCL21"]





sc.pl.dotplot(adata,
              marker_genes,
              groupby='cell_type_am', use_raw=True,standard_scale="var")





sc.pl.umap(adata,color="cell_type_am",legend_loc="on data", title='', frameon=False,legend_fontsize="x-small")





#reverse adata.X raw
adata.X=adata.layers["raw_counts"].copy()





adata.write(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated.h5ad")


# # Save separate objects for each cell_type for downstream clustering




for i in set(adata.obs.cell_type_am):
    tmp=adata[adata.obs.cell_type_am.isin([i]),].copy()
    tmp.write(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_"+i+".h5ad")
    print("Saved "+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_"+i+".h5ad")
    del tmp

