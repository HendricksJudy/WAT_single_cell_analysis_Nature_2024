#!/usr/bin/env python
# coding: utf-8




import pandas as pd
import numpy as np
import gc
import scanpy as sc
import loompy as lp
import re
import matplotlib.pyplot as plt

import scipy
import scipy.sparse
import seaborn as sns

from sklearn.preprocessing import StandardScaler





sc.settings.set_figure_params(dpi=300,fontsize=10,dpi_save=300)





scaler = StandardScaler()





def name_of_object(arg):
    # check __name__ attribute (functions)
    try:
        return arg.__name__
    except AttributeError:
        pass

    for name, value in globals().items():
        if value is arg and not name.startswith('_'):
            return name


# # Load data




path="path_in/"
INPUT = path+"swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad"





adata = sc.read(INPUT)





adata





adata.X=adata.layers["log1p_counts"]





del adata.layers
gc.collect()





adata=adata[adata.obs.dataset.isin(["scott"])]


# # Stress score




stress_genes=pd.read_csv("stress_state_markers_vs_basal.csv",index_col=0)





stress_genes=pd.DataFrame(stress_genes[(stress_genes.lfc>=0.1)&(~stress_genes.cell_state.isin(["Mu5"]))&(stress_genes.adj_pval<=0.05)].Gene.value_counts()) # set threshold and exlude Mu5 due to low power





stress_genes=np.asarray(stress_genes[stress_genes.Gene==4].index) # keep genes that appear in all stressed cell states





sc.tl.score_genes(adata, stress_genes, ctrl_size=50, score_name='stress_score', use_raw=None)





adata.obs["ct_condition"]=adata.obs["cell_type_am_fine"].astype(str)+"_"+adata.obs["condition2"].astype(str)





# Normalize the 'score' variable in adata.obs
score = adata.obs['stress_score']

# Perform z-score normalization (standardize the 'score' variable)
normalized_score = (score - score.mean()) / score.std()

adata.obs['normalized_score'] = normalized_score





# for cell type
# Extract the relevant categorical and continuous variables
df = adata.obs[['condition2', 'cell_type_am_fine', 'normalized_score']].copy()

# Create the pivot table: rows = 'condition', columns = 'cluster', values = normalized_score
pivot_data = df.pivot_table(values='normalized_score', index='condition2', columns='cell_type_am_fine', aggfunc='mean')





del pivot_data["Unassigned"]





# for cell state
# Extract the relevant categorical and continuous variables
df = adata.obs[['condition2', 'cell_state_am', 'normalized_score']].copy()

# Create the pivot table: rows = 'condition', columns = 'cluster', values = normalized_score
pivot_data = df.pivot_table(values='normalized_score', index='condition2', columns='cell_state_am', aggfunc='mean')





del pivot_data["Unassigned"]


# # Subset to adipocytes for next section (can be done for any cell)




adata=adata[adata.obs.cell_type_am_fine.isin(["Adipocytes"])]





gc.collect()


# # Senescence scores

# ### SenMayo




senMayo = pd.read_csv("senescence_SenMayo_human.csv") 





senMayo = [x for x in senMayo["Gene(human)"] if x in adata.var_names]





sc.tl.score_genes(adata, senMayo, ctrl_size=50, score_name='senMayo_score', use_raw=None)





gc.collect()


# ### Fridman Senescence list




sen=pd.read_csv("geneset_FRIDMAN_SENESCENCE_UP_HUMAN.csv")





sen_list=[x for x in sen.iloc[:, 0] if x in np.asarray(adata.var.index)]





sc.tl.score_genes(adata, sen_list, ctrl_size=50, score_name='Fridman_score', use_raw=None)





gc.collect()


# ### SASP list




sen=pd.read_csv("SASP_Signatures_list.csv",sep="\t",header=None)





sen_list=[x for x in sen.iloc[:, 0] if x in np.asarray(adata.var.index)]





sc.tl.score_genes(adata, sen_list, ctrl_size=50, score_name='SASP_score', use_raw=None)





gc.collect()


# ### OIS list




sen=pd.read_csv("OIS_Gene_List.csv",header=None)





sen_list=[x for x in sen.iloc[:, 0] if x in np.asarray(adata.var.index)]





sc.tl.score_genes(adata, sen_list, ctrl_size=50, score_name='OIS_score', use_raw=None)





gc.collect()


# ### Broad list




sen=pd.read_csv("Broad.MSigDB.AllSenescent.csv")





for j in sen.columns:
    sen_list=[x for x in sen[j] if x in np.asarray(adata.var.index)]
    score_name="Broad_"+j
    sc.tl.score_genes(adata, sen_list, ctrl_size=50, score_name=score_name, use_raw=None)
    gc.collect()


# ### All sen scores




sen_scores=['senMayo_score', 'Fridman_score', 'SASP_score', 'OIS_score',
       'Broad_BIOCARTA_TEL_PATHWAY', 'Broad_CHICAS_RB1_TARGETS_CONFLUENT',
       'Broad_CHICAS_RB1_TARGETS_GROWING',
       'Broad_CHICAS_RB1_TARGETS_LOW_SERUM',
       'Broad_CHICAS_RB1_TARGETS_SENESCENT',
       'Broad_COURTOIS_SENESCENCE_TRIGGERS', 'Broad_DEMAGALHAES_AGING_DN',
       'Broad_DEMAGALHAES_AGING_UP', 'Broad_FRIDMAN_IMMORTALIZATION_DN',
       'Broad_FRIDMAN_SENESCENCE_DN', 'Broad_FRIDMAN_SENESCENCE_UP',
       'Broad_GO_CELL_AGING', 'Broad_GO_MULTICELLULAR_ORGANISM_AGING',
       'Broad_GO_STRESS_INDUCED_PREMATURE_SENESCENCE',
       'Broad_KAMMINGA_EZH2_TARGETS', 'Broad_KAMMINGA_SENESCENCE',
       'Broad_KEGG_P53_SIGNALING_PATHWAY',
       'Broad_KUMAMOTO_RESPONSE_TO_NUTLIN_3A_DN',
       'Broad_KUMAMOTO_RESPONSE_TO_NUTLIN_3A_UP',
       'Broad_LINDVALL_IMMORTALIZED_BY_TERT_DN',
       'Broad_LINDVALL_IMMORTALIZED_BY_TERT_UP',
       'Broad_ONGUSAHA_BRCA1_TARGETS_DN',
       'Broad_ONGUSAHA_BRCA1_TARGETS_UP', 'Broad_ONGUSAHA_TP53_TARGETS',
       'Broad_REACTOME_CELLULAR_SENESCENCE',
       'Broad_REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE',
       'Broad_REACTOME_ONCOGENE_INDUCED_SENESCENCE',
       'Broad_REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE',
       'Broad_ROETH_TERT_TARGETS_DN', 'Broad_ROETH_TERT_TARGETS_UP',
       'Broad_TANG_SENESCENCE_TP53_TARGETS_DN',
       'Broad_TANG_SENESCENCE_TP53_TARGETS_UP',
       'Broad_VARELA_ZMPSTE24_TARGETS_DN',
       'Broad_VARELA_ZMPSTE24_TARGETS_UP',
       'Broad_WP_SENESCENCE_AND_AUTOPHAGY_IN_CANCER']


# # Metabolic scores




pik3=adata.var.index[adata.var.index.str.startswith("PIK3")]





tca = np.unique(["ACO2", "CS", "FH", "MDH1", "OGDH", "PDHA1", "SDHC", "SUCLG1","SDHA","SDHB","SDHD","SUCLA2","IDH1","IDH2","OGDH","OGDHL","MDH1B","MDH2"])
dnl = np.unique(["MLXIPL", "NR1H3","NR1H2", "SREBF1", "FASN", "ACACA","ACLY","ELOVL6","SCD","ACACB","MECR","OXSM"])
insulin = np.unique(np.append(["INSR","IRS1","IRS2","PDE3B","MTOR","SLC2A4","SLC2A1","AKT1","AKT2","AKT3"],pik3))
lipolysis = np.unique(["PNPLA2","LIPE","MGLL","CIDEA","CAVIN1","CAV1"])
betaox = np.unique(["ACSS2","ACADM","ECHS1","HADHA","HADHB","ACAT1","ACAA1","ACAA2","ECH1","ECHS1","ECI1","ECI2","EHHADH","ACACB"])
traffick = np.unique(["CD36","LPL","FABP4"])
elongation = np.unique(["ELOVL1","ELOVL2","ELOVL3","ELOVL4","ELOVL5","ELOVL6","ELOVL7","HACD1","HACD2","HACD3","HACD4","HIGD2A","HSD17B12","TECR","TECRL","OXSM","PECR"])
esterification = np.unique(["GPAM", "AGPAT1", "AGPAT2", "AGPAT3", "AGPAT4", "AGPAT5", "DGAT1", "DGAT2", "PLPP1", "PLPP2", "PLPP3", "PLPP4" ,"PLPP5","ACAT1","ACAT2","GPAM","GPAT2","GPAT3","ACSS2","MOGAT1","MOGAT2","MOGAT3"])
lipid_droplet = np.unique(["PLIN1","PLIN2","PLIN3","PLIN4","PLIN5","CIDEA","CIDEC"])



creat_thermo = np.unique(["SLC6A8","GATM","GAMT","CKMT1A","CKMT1B","CKMT2","CKB"])
cal_therm = np.unique(["ATP2A2","RYR2","ITPR1","ITPR2","ITPR3","ATP2A1"])
adap_thermo = np.unique(["EBF2","ESRRG","PPARGC1A","PPARGC1B","PRKAA1","PRKAG2","PERM1","ATP2A1","ATP2A2","ATP2A3","ESRRA","DIO2","PRDM16"])


bcaa = np.unique(['ABAT', 'ACAA2', 'ACADM', 'IVD', 'ACAD8', 'ECHS1', 'EHHADH',
       'HIBADH', 'BCAT1', 'BCAT2', 'MCCC1', 'AUH', 'MLYCD', 'MCEE', 'MUT',
       'ALDH6A1', 'OXCT1', 'BCKDHA', 'PCCA', 'ACADS', 'DBT', 'ALDH1B1',
       'BCKDHB', 'ACAT1', 'ACADSB', 'HADHA', 'HADH', 'MCCC2', 'OXCT2B',
       'PCCB', 'ALDH2', 'HADHB', 'HSD17B10', 'ALDH3A2', 'DLD', 'ALDH7A1',
       'ALDH9A1'])





scores=[tca,dnl,insulin,lipolysis,betaox,traffick,elongation,esterification,lipid_droplet,creat_thermo,cal_therm,adap_thermo,bcaa]





for k in scores:
    name=name_of_object(k)
    sc.tl.score_genes(adata, k, ctrl_size=50, score_name=name, use_raw=None)





housekeeping=["RRN18S", "ACTB", "GAPDH", "PGK1", "PPIA", "RPL13A", "RPLP0", "ARBP", "B2M", "YWHAZ", "SDHA", "TFRC", "GUSB", "HMBS", "HPRT1", "TBP",
             "EEF2","MRFAP1","RAB7A","RAB1B","HNRNPA2B1","PFDN5","GABARAP","RAB11B","CSNK2B","RHOA","BSG","PCBP1","UBC","HSP90AB1","AP2M1","EDF1","AES","MLF2","PSAP","CD59"]





housekeeping=[x for x in housekeeping if x in adata.var.index]





sc.tl.score_genes(adata, housekeeping, ctrl_size=50, score_name='Housekeeping', use_raw=None)





scores=["tca","dnl","insulin","lipolysis","betaox","traffick","elongation","esterification","lipid_droplet","creat_thermo","cal_therm","adap_thermo","bcaa","Housekeeping","stress_score"]


# # Correlate senescence and metabolic scores




selected_sen_scores=["Broad_FRIDMAN_SENESCENCE_UP","Broad_GO_CELL_AGING",
                "Broad_KEGG_P53_SIGNALING_PATHWAY","Broad_REACTOME_CELLULAR_SENESCENCE","Broad_WP_SENESCENCE_AND_AUTOPHAGY_IN_CANCER"]





all_scores=np.append(scores,selected_sen_scores)





df=adata.obs[all_scores]





df





corr_matrix=df.corr()


# ### delta correlation




obese_df=df[adata.obs.condition2.isin(["Obese"])]
lean_df=df[adata.obs.condition2.isin(["Lean"])]
wl_df=df[adata.obs.condition2.isin(["Weightloss"])]





obese_lean_delta=obese_df-lean_df.mean()
obese_lean_delta=obese_lean_delta.corr()





obese_df["donor_id"]=adata[adata.obs.condition2.isin(["Obese"])].obs.donor_id
wl_df["donor_id"]=adata[adata.obs.condition2.isin(["Weightloss"])].obs.donor_id





# paired delta for OB-WL
final_df = None  # Initialize final_df to None before the loop
for i in obese_df.donor_id.unique():
    # Select data for the current donor in both obese and wl dataframes
    tmp_df_ob = obese_df[obese_df.donor_id == i].drop(columns=['donor_id'])  # Drop 'donor_id' for the subtraction
    tmp_df_wl = wl_df[wl_df.donor_id == i].drop(columns=['donor_id'])  # Same here
    
    # Subtract the mean of wl_df for each donor from the obese dataframe
    tmp_df_ob = tmp_df_ob - tmp_df_wl.mean()

    # Concatenate the results into final_df
    if final_df is None:
        final_df = tmp_df_ob.copy()  # Initialize final_df with the first donor's data
    else:
        final_df = pd.concat([final_df, tmp_df_ob], ignore_index=True)  # Concatenate for subsequent donors
final_df=final_df.corr()







