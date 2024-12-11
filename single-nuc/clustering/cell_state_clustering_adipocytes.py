#!/usr/bin/env python
# coding: utf-8




#import all tools  - if not installed please install them using "pip install toolname --user" or following the instructions on their GitHub page.
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import matplotlib.pyplot as plt
import os
import bbknn as bk
import anndata 
import harmonypy as hm
get_ipython().run_line_magic('matplotlib', 'inline')



sc.settings.set_figure_params(dpi=300,fontsize=10,dpi_save=300)



# Define required functions

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
        #incluster_rows, incluster_cols=adata.X.nonzero()
        #https://stackoverflow.com/questions/3797158/counting-non-zero-elements-within-each-row-and-within-each-column-of-a-2d-numpy
        #incluster_counts = [len(np.where(incluster_cols==i)[0]) for i in range(0,adata.shape[1])]
        incluster_counts=np.diff(adata.X.tocsc().indptr)
        percent_expression_incluster=np.array(incluster_counts)/sizeof_adata 
        fraction_in_cluster_matrix[cluster_number[0]]=percent_expression_incluster
        #outcluster_rows, outcluster_cols=adata_outcluster.X.nonzero()
        #outcluster_counts = [len(np.where(outcluster_cols==i)[0]) for i in range(0,adata_outcluster.shape[1])]
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


# # Load Dataset




path_in="path_in/"





seed=5





adata=sc.read_h5ad(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_Adipocytes.h5ad")





adata






del adata.uns["log1p"]
adata.layers["raw_counts_state"]=adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
adata.layers["norm_counts_state"]=adata.X.copy()
sc.pp.log1p(adata)
adata.layers["log1p_counts_state"]=adata.X.copy()
adata.raw=adata.copy()





#verify if adata.X is log transformed
adata.X.max()
#it is log transford





#verify if adata.raw.X is log transformed
adata.raw.X.max()
# it is log transformed


# # Initial clustering to indentify "doublet" cell states




sc.pp.regress_out(adata,["mt.percent","nCount_RNA","ribo.percent"])





sc.pp.highly_variable_genes(adata)





sc.tl.pca(adata, svd_solver='arpack',n_comps=40,random_state=seed,use_highly_variable=True)





pca = adata.obsm['X_pca']
batch = adata.obs['sample']
meta_data = adata.obs

ho = hm.run_harmony(pca, meta_data, ['sample'], epsilon_harmony = -float('Inf'), max_iter_kmeans = 15, tau=5, max_iter_harmony = 25)





res = pd.DataFrame(ho.Z_corr)
res=res.T
adata.obsm['X_pca'] = res.values





bk.bbknn(adata,n_pcs=40,local_connectivity=1,batch_key="sample",neighbors_within_batch=1)





sc.tl.umap(adata,spread=0.8,min_dist=0.15)





# low resolution to identify only major cell types
sc.tl.leiden(adata, key_added="cell_state_am",resolution=1.2,random_state=seed,n_iterations=-1)





sc.pl.umap(adata,color="cell_state_am")





sc.pl.umap(adata,color="sample")





sc.pl.umap(adata,color=["dataset","condition"])





sc.pl.umap(adata,color=["nCount_RNA","mt.percent","ribo.percent"])





marker_genes=["ADIPOQ","GPAM","DCN","PDGFRA","PTPRC","MRC1","CD68","VWF","CDH5","CD247","NKG7","MYH11","KCNJ8","KIT","CPA3"]





sc.pl.dotplot(adata,
              marker_genes,
              groupby='cell_state_am', use_raw=True,standard_scale="var")





sc.pl.violin(adata,["scrub_score"],groupby="cell_state_am")





# Cluster 7 likely countain doublets. They have a higher scrublet scores and express pan-markers of other lineages.
# Call all this unassigned to facilitate removal and later identification
#rename clusters to respective cell type
adata.obs.cell_state_am=adata.obs.cell_state_am.replace(
    {"7":"Unassigned"
    })





#revert adata.X to non regressed state
adata.X=adata.raw.X





adata.write(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_Adipocytes.h5ad")


# # Remove doublets and re-cluster




adata_clean=adata[~adata.obs.cell_state_am.isin(["Unassigned"])].copy()





adata_clean





sc.pp.regress_out(adata_clean,["mt.percent","nCount_RNA","ribo.percent"])





sc.pp.highly_variable_genes(adata_clean)





sc.tl.pca(adata_clean, svd_solver='arpack',n_comps=40,random_state=seed,use_highly_variable=True)





pca = adata_clean.obsm['X_pca']
batch = adata_clean.obs['sample']
meta_data = adata_clean.obs

ho = hm.run_harmony(pca, meta_data, ['sample'], epsilon_harmony = -float('Inf'), max_iter_kmeans = 15, tau=5, max_iter_harmony = 25)





res = pd.DataFrame(ho.Z_corr)
res=res.T
adata_clean.obsm['X_pca'] = res.values





bk.bbknn(adata_clean,n_pcs=40,local_connectivity=1,batch_key="sample",neighbors_within_batch=1)





sc.tl.umap(adata_clean,spread=0.8,min_dist=0.15)





# low resolution to identify only major cell types
sc.tl.leiden(adata_clean, key_added="cell_state_am",resolution=0.95,random_state=seed,n_iterations=-1)





sc.pl.umap(adata_clean,color="cell_state_am")





sc.pl.umap(adata_clean,color="sample")





sc.pl.umap(adata_clean,color=["dataset","condition"])





sc.pl.umap(adata_clean,color=["nCount_RNA","mt.percent","ribo.percent"])





marker_genes=["ADIPOQ","GPAM","DCN","PDGFRA","PTPRC","MRC1","CD68","VWF","CDH5","CD247","NKG7","MYH11","KCNJ8","KIT","CPA3"]





sc.pl.dotplot(adata_clean,
              marker_genes,
              groupby='cell_state_am', use_raw=True,standard_scale="var")





sc.pl.violin(adata_clean,["scrub_score"],groupby="cell_state_am")





sc.tl.rank_genes_groups(adata_clean, "cell_state_am", method="t-test_overestim_var",use_raw=True)
sc.tl.dendrogram(adata_clean,groupby="cell_state_am")





sc.pl.rank_genes_groups_dotplot(adata_clean, n_genes=5, groupby="cell_state_am",mean_only_expressed=True,standard_scale="var")





#this adds the expression in and expression out information to adata and creates a table that can be exported

expression=filter_by_expression(adata_clean, 'cell_state_am', key='rank_genes_groups')





expression





#Dotplot with filtering options for unbiased markers
#Arguments for out cluster are <=,i.e, a percent_out=1 means no filtering will occur based on percentage outside the cluster
#Arguments for in cluster are >=,i.e, a percent_in=0 means no filtering will occur based on percentage inside the cluster


sc.pl.dotplot(adata_clean,
              makevector_topmarkers(adata_clean, 
                                    5, ldFC_cutoff=0.5, pvals_adj_cutoff=10**-3, expression_in=0.1,expression_out=1,percentage_out=1,percentage_in=0.3),
              groupby='cell_state_am', use_raw=True,standard_scale="var")





marker_genes={
    "Emont":["GALNT13","TNFSF10","PNPLA3","GRIA4","PGAP1","EBF2","AGMO"],
    "Stress":["ATF3","FOS"]
}





sc.pl.dotplot(adata_clean,
              marker_genes,
              groupby='cell_state_am', use_raw=True,standard_scale="var")





#Merge clusters?
#clusters 1 and 3 have similar signatures. Are there differences:
subset=adata_clean[adata_clean.obs.cell_state_am.isin(["1","3"])].copy()
sc.tl.rank_genes_groups(subset, "cell_state_am", method="t-test_overestim_var",use_raw=True)
sc.tl.dendrogram(subset,groupby="cell_state_am")
sc.pl.rank_genes_groups_dotplot(subset, n_genes=20, groupby="cell_state_am",mean_only_expressed=True)
#No extensive differences, will merge the clusters





#merge cluster 2 and 4
adata_clean.obs.cell_state_am=adata_clean.obs.cell_state_am.replace(
    {"3":"1"
    })





sc.pl.umap(adata_clean,color="cell_state_am")





sc.tl.rank_genes_groups(adata_clean, "cell_state_am", method="t-test_overestim_var",use_raw=True)
sc.tl.dendrogram(adata_clean,groupby="cell_state_am")





sc.pl.rank_genes_groups_dotplot(adata_clean, n_genes=5, groupby="cell_state_am",mean_only_expressed=True,standard_scale="var")





#this adds the expression in and expression out information to adata and creates a table that can be exported

expression=filter_by_expression(adata_clean, 'cell_state_am', key='rank_genes_groups')





expression





#Dotplot with filtering options for unbiased markers
#Arguments for out cluster are <=,i.e, a percent_out=1 means no filtering will occur based on percentage outside the cluster
#Arguments for in cluster are >=,i.e, a percent_in=0 means no filtering will occur based on percentage inside the cluster


sc.pl.dotplot(adata_clean,
              makevector_topmarkers(adata_clean, 
                                    5, ldFC_cutoff=0.5, pvals_adj_cutoff=10**-3, expression_in=0.1,expression_out=1,percentage_out=1,percentage_in=0.3),
              groupby='cell_state_am', use_raw=True,standard_scale="var")





marker_genes={
    "Emont":["GALNT13","TNFSF10","PNPLA3","GRIA4","PGAP1","EBF2","AGMO"],
    "Stress":["ATF3","FOS"]
}





sc.pl.dotplot(adata_clean,
              marker_genes,
              groupby='cell_state_am', use_raw=True,standard_scale="var")





#Merge clusters?
#clusters 2 and 5 have no clear markers and could represent a basal state. Are there differences:
subset=adata_clean[adata_clean.obs.cell_state_am.isin(["2","5"])].copy()
sc.tl.rank_genes_groups(subset, "cell_state_am", method="t-test_overestim_var",use_raw=True)
sc.tl.dendrogram(subset,groupby="cell_state_am")
sc.pl.rank_genes_groups_dotplot(subset, n_genes=20, groupby="cell_state_am",mean_only_expressed=True)
#No extensive differences, will merge the clusters





#merge cluster 2 and 4
adata_clean.obs.cell_state_am=adata_clean.obs.cell_state_am.replace(
    {"5":"2"
    })





sc.pl.umap(adata_clean,color="cell_state_am")





sc.tl.rank_genes_groups(adata_clean, "cell_state_am", method="t-test_overestim_var",use_raw=True)
sc.tl.dendrogram(adata_clean,groupby="cell_state_am")





sc.pl.rank_genes_groups_dotplot(adata_clean, n_genes=5, groupby="cell_state_am",mean_only_expressed=True,standard_scale="var")





#this adds the expression in and expression out information to adata and creates a table that can be exported

expression=filter_by_expression(adata_clean, 'cell_state_am', key='rank_genes_groups')





expression





#Dotplot with filtering options for unbiased markers
#Arguments for out cluster are <=,i.e, a percent_out=1 means no filtering will occur based on percentage outside the cluster
#Arguments for in cluster are >=,i.e, a percent_in=0 means no filtering will occur based on percentage inside the cluster


sc.pl.dotplot(adata_clean,
              makevector_topmarkers(adata_clean, 
                                    5, ldFC_cutoff=0.5, pvals_adj_cutoff=10**-3, expression_in=0.1,expression_out=1,percentage_out=1,percentage_in=0.3),
              groupby='cell_state_am', use_raw=True,standard_scale="var")





marker_genes={
    "Emont":["GALNT13","TNFSF10","PNPLA3","GRIA4","PGAP1","EBF2","AGMO"],
    "Stress":["ATF3","FOS"]
}





sc.pl.dotplot(adata_clean,
              marker_genes,
              groupby='cell_state_am', use_raw=True,standard_scale="var")





adata_clean.obs["cell_state_am_long"]=adata_clean.obs["cell_state_am"]
adata_clean.obs.cell_state_am_long=adata_clean.obs.cell_state_am_long.replace(
    {"2":"AD1_basal",
     "0":"AD2_GRIA4_high",
     "1":"AD3_stressed",
     "4":"AD4_intermediate",
     "6":"AD5_PNPLA3hi",
     "7":"AD6_SLC14A2hi",
     "8":"AD7_PGAP1hi",
     "9":"AD8_AGMOhi"
    })





sc.pl.umap(adata_clean,color="cell_state_am_long")





#rename numbers on cell state column to a simple designation
adata_clean.obs.cell_state_am=adata_clean.obs.cell_state_am.replace(
    {"2":"AD1",
     "0":"AD2",
     "1":"AD3",
     "4":"AD4",
     "6":"AD5",
     "7":"AD6",
     "8":"AD7",
     "9":"AD8"
    })





sc.pl.umap(adata_clean,color="cell_state_am")





#revert adata.X to non regressed state
adata_clean.X=adata_clean.raw.X





# save object with no doublets
adata_clean.write(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_Adipocytes_noDoublets.h5ad")





sc.tl.rank_genes_groups(adata_clean, "cell_state_am_long", method="t-test_overestim_var",use_raw=True)
sc.tl.dendrogram(adata_clean,groupby="cell_state_am_long")





sc.pl.rank_genes_groups_dotplot(adata_clean, n_genes=5, groupby="cell_state_am_long",mean_only_expressed=True,standard_scale="var")


# # Add cell state information to dataset with doublets




adata=sc.read_h5ad(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_Adipocytes.h5ad")
adata_clean=sc.read_h5ad(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_Adipocytes_noDoublets.h5ad")





df=pd.concat([adata_clean.obs.cell_state_am.rename("simple").astype(str),adata_clean.obs.cell_state_am_long.rename("long").astype(str),adata.obs.cell_state_am.rename("doublet").astype(str)],axis=1,names=["clean","doublet"])





df





df['simple'][df["doublet"].isin(['Unassigned'])] = "Unassigned"
df['long'][df["doublet"].isin(['Unassigned'])] = "Unassigned"





df





adata.obs.cell_state_am=df["simple"]





adata.obs["cell_state_am_long"]=df["long"]





sc.pl.umap(adata,color="cell_state_am_long")


# ### Add UMAP coordinates to adata.obs




coord=pd.DataFrame(adata.obsm["X_umap"],columns=["UMAP1_doublets","UMAP2_doublets"],index=adata.obs.index)





coord





coord2=pd.DataFrame(adata_clean.obsm["X_umap"],columns=["UMAP1_clean","UMAP2_clean"],index=adata_clean.obs.index)





coord2





coord=pd.concat([coord,coord2],axis=1)





coord





adata.obs["UMAP1_doublets"]=coord["UMAP1_doublets"]
adata.obs["UMAP2_doublets"]=coord["UMAP2_doublets"]
adata.obs["UMAP1_clean"]=coord["UMAP1_clean"]
adata.obs["UMAP2_clean"]=coord["UMAP2_clean"]





adata





#Write object that includes doublets
adata.write(path_in+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_Adipocytes.h5ad")







