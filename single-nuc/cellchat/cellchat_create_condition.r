
library(reticulate)
library(CellChat)
library(patchwork)
library(stringr)
options(stringsAsFactors = FALSE)


ad <- import("anndata", convert = FALSE)

path_in="/cellchat/input_files/"

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

ad_lean <- ad$read_h5ad(paste0(path_in,"global_Lean_for_cellchat_no_lymphatic_bcells_kit_scott_only_stringent_subsampled.h5ad"))

data=as.matrix(py_to_r(ad_lean$X))

gc()

data.input <- t(data)

rownames(data.input) <- rownames(py_to_r(ad_lean$var))
colnames(data.input) <- rownames(py_to_r(ad_lean$obs))

meta.data <- py_to_r(ad_lean$obs)
meta <- meta.data

#cellchat labels cannot have "_". Best avoid all other special characters to.
meta$cell_type_am_fine <- str_replace(meta$cell_type_am_fine, "_", ".")
meta$cell_type_am_fine <- str_replace(meta$cell_type_am_fine, "-", ".")
meta$cell_type_am_fine <- str_replace(meta$cell_type_am_fine, "\\+", "")

meta$cell_type_am_fine=factor(meta$cell_type_am_fine)

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type_am_fine")

# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file =paste0(path_out,"cellchat_lean_no_lymphatic_bcells_kit_scott_only_ct_fine_stringent_subsampled.rds"))


