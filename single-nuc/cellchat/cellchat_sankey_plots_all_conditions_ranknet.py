#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import numpy as np
import gc
import scanpy as sc
import re
import colorcet as cc
import seaborn as sns
import matplotlib.colors
from PIL import ImageColor
import plotly.graph_objects as go
import plotly.offline as pyo
pyo.init_notebook_mode()
import plotly.io as pio
pio.renderers.default = 'iframe'





path_cellchat="/cellchat/output_files/"





ranknet=pd.read_csv(path_cellchat+"rankNet_vs_obese_ct_fine_pairwise_scott_only_subsampled_stringent_degs_no_lymphatic_bcells_kit.csv",index_col=0) #obtained from cellchat_multiple_comparisons.r





ranknet





ranknet.contrib.max()


# # Weightloss

# # Nodes




df=ranknet[(ranknet.Condition.isin(["Weightloss"]))&(~ranknet["contrib"].isin([np.nan]))&(ranknet["pval"]<0.05)]





df["source_2"]=df["Sender"].astype(str)+"_s"
df["target_2"]=df["Receiver"].astype(str)+"_t"





celltypes=np.append(df["target_2"].unique(),df["source_2"].unique())
celltypes=np.unique(celltypes)
pathways=df["Pathway"].unique()
nodes=np.append(celltypes,pathways)
nodes=pd.DataFrame(nodes)





nodes.columns=["Label"]
nodes["ID"]=nodes.index





nodes[["Label_real","type"]]=nodes["Label"].str.split("_",expand=True)





nodes_unique=nodes.drop_duplicates("Label_real")





nodes_unique["color"]=cc.glasbey[:nodes_unique.shape[0]]





color_dict=pd.Series(nodes_unique["color"].values,index=nodes_unique["Label_real"]).to_dict()





nodes["color"]=nodes["Label_real"].replace(color_dict)





rgb_colors=[ImageColor.getcolor(x, "RGBA") for x in nodes["color"]]





nodes["rgb_color"]=rgb_colors





nodes["rgb_color"] = "rgba"+ nodes["rgb_color"].astype(str)
nodes["rgb_color_link"] =  nodes["rgb_color"].apply(lambda x: x.replace('255)','0.5)'))
                                                    





replace_dict=pd.Series(nodes["ID"].values,index=nodes["Label"]).to_dict()





color_dict=pd.Series(nodes["rgb_color_link"].values,index=nodes["ID"]).to_dict()


# # Links




df["source_pathway"]=df["source_2"].astype(str)+"__"+df["Pathway"].astype(str)





source=pd.DataFrame(df["contrib"].groupby(df["source_pathway"]).mean())





source=pd.concat([source,df["source_pathway"].value_counts()],axis=1)





source.columns=["contrib","Value"]





df["pathway_target"]=df["Pathway"].astype(str)+"__"+df["target_2"].astype(str)





target=pd.DataFrame(df["contrib"].groupby(df["pathway_target"]).mean())





target=pd.concat([target,df["pathway_target"].value_counts()],axis=1)





target.columns=["contrib","Value"]





links=pd.concat([source,target],axis=0)





links["pair"]=links.index





links[["source","target"]]=links["pair"].str.split("__",expand=True)
del links["pair"]





links.sort_values(by="Value",ascending=False).head()





links["source"]=links["source"].replace(replace_dict)
links["target"]=links["target"].replace(replace_dict)





links.sort_values(by="Value",ascending=False).head()





links["link_color"]=links["source"].replace(color_dict)





links["abs_value"]=abs(links["contrib"])





orange="#EF5600"
darkred="#d60012"
blue="#5689EF"
darkblue="#0017eb"
gray="#a5a6ad"





links["color_condition"]=orange
links["color_condition"]=np.where((links["contrib"]>0),blue,links["color_condition"])





links["alpha_scale"]=links["abs_value"]/links["abs_value"].max()





rgb_colors=[ImageColor.getcolor(x, "RGBA") for x in links["color_condition"]]





links["color_condition"]=rgb_colors





links["color_condition"] = "rgba"+ links["color_condition"].astype(str)
links["color_condition"] =  links["color_condition"].apply(lambda x: x.replace('255)',''))
links["color_condition"]=links["color_condition"].astype(str)+links["alpha_scale"].astype(str)+")"


# # Plot




fig = go.Figure(data=[go.Sankey(
    valueformat = ".0f",
    # Define nodes
    node = dict(
      pad = 15,
      thickness = 15,
      line = dict(color = "black", width = 0.5),
      label =  nodes['Label_real'],
      color =  nodes['rgb_color']
    ),
    # Add links
    link = dict(
      source =  links['source'],
      target =  links['target'],
      value =  links['Value'],
      color =  links['color_condition']
))])

fig.update_layout(
    autosize=False,
    width=800,
    height=1200,
    font_family="arial"
)
fig.write_image("sankey_weightloss_ranknet.pdf")
fig.show()


# # Lean

# # Nodes




df=ranknet[(ranknet.Condition.isin(["Lean"]))&(~ranknet["contrib"].isin([np.nan]))&(ranknet["pval"]<0.05)]





df["source_2"]=df["Sender"].astype(str)+"_s"
df["target_2"]=df["Receiver"].astype(str)+"_t"





celltypes=np.append(df["target_2"].unique(),df["source_2"].unique())
celltypes=np.unique(celltypes)
pathways=df["Pathway"].unique()
nodes=np.append(celltypes,pathways)
nodes=pd.DataFrame(nodes)





nodes.columns=["Label"]
nodes["ID"]=nodes.index





nodes[["Label_real","type"]]=nodes["Label"].str.split("_",expand=True)





nodes_unique=nodes.drop_duplicates("Label_real")





nodes_unique["color"]=cc.glasbey[:nodes_unique.shape[0]]





color_dict=pd.Series(nodes_unique["color"].values,index=nodes_unique["Label_real"]).to_dict()





nodes["color"]=nodes["Label_real"].replace(color_dict)





rgb_colors=[ImageColor.getcolor(x, "RGBA") for x in nodes["color"]]





nodes["rgb_color"]=rgb_colors





nodes["rgb_color"] = "rgba"+ nodes["rgb_color"].astype(str)
nodes["rgb_color_link"] =  nodes["rgb_color"].apply(lambda x: x.replace('255)','0.5)'))
                                                    





replace_dict=pd.Series(nodes["ID"].values,index=nodes["Label"]).to_dict()





color_dict=pd.Series(nodes["rgb_color_link"].values,index=nodes["ID"]).to_dict()


# # Links




df["source_pathway"]=df["source_2"].astype(str)+"__"+df["Pathway"].astype(str)





source=pd.DataFrame(df["contrib"].groupby(df["source_pathway"]).mean())





source=pd.concat([source,df["source_pathway"].value_counts()],axis=1)





source.columns=["contrib","Value"]





df["pathway_target"]=df["Pathway"].astype(str)+"__"+df["target_2"].astype(str)





target=pd.DataFrame(df["contrib"].groupby(df["pathway_target"]).mean())





target=pd.concat([target,df["pathway_target"].value_counts()],axis=1)





target.columns=["contrib","Value"]





links=pd.concat([source,target],axis=0)





links["pair"]=links.index





links[["source","target"]]=links["pair"].str.split("__",expand=True)
del links["pair"]





links.sort_values(by="Value",ascending=False).head()





links["source"]=links["source"].replace(replace_dict)
links["target"]=links["target"].replace(replace_dict)





links.sort_values(by="Value",ascending=False).head()





links["link_color"]=links["source"].replace(color_dict)





links["abs_value"]=abs(links["contrib"])





orange="#EF5600"
darkred="#d60012"
blue="#5689EF"
darkblue="#0017eb"
gray="#a5a6ad"





links["color_condition"]=orange
links["color_condition"]=np.where((links["contrib"]>0),blue,links["color_condition"])





links["alpha_scale"]=links["abs_value"]/links["abs_value"].max()





rgb_colors=[ImageColor.getcolor(x, "RGBA") for x in links["color_condition"]]





links["color_condition"]=rgb_colors





links["color_condition"] = "rgba"+ links["color_condition"].astype(str)
links["color_condition"] =  links["color_condition"].apply(lambda x: x.replace('255)',''))
links["color_condition"]=links["color_condition"].astype(str)+links["alpha_scale"].astype(str)+")"


# # Plot




fig = go.Figure(data=[go.Sankey(
    valueformat = ".0f",
    # Define nodes
    node = dict(
      pad = 15,
      thickness = 15,
      line = dict(color = "black", width = 0.5),
      label =  nodes['Label_real'],
      color =  nodes['rgb_color']
    ),
    # Add links
    link = dict(
      source =  links['source'],
      target =  links['target'],
      value =  links['Value'],
      color =  links['color_condition']
))])

fig.update_layout(
    autosize=False,
    width=800,
    height=1200,
    font_family="arial"
)
fig.write_image("sankey_lean_ranknet.pdf")
fig.show()







