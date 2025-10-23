# Figure 2 Single-Nucleus Analysis Pipeline

This document describes the end-to-end workflow for generating the panels shown in Figure 2a, 2e, and 2g of Miranda *et al.* (Nature 2024). The steps below assume access to the integrated single-nucleus dataset `swat_all_seurat_integration_scott_ref_stringent_rna.h5ad` produced by the global integration pipeline. Paths mentioned below are relative to the repository root unless otherwise noted.

> **New:** the automated runner in [`single-nuc/pipelines/figure2`](./figure2) wraps the steps below so the full set of Figure 2 panels can be regenerated with a single command.

```bash
python single-nuc/pipelines/figure2/cli.py config/figure2.yaml
```

An example configuration can be used as a template:

```yaml
integrated_h5ad: data/swat_all_seurat_integration_scott_ref_stringent_rna.h5ad
output_dir: results/figure2
preprocess:
  cluster_cell_type_map:
    "0": Adipocyte
    "1": Myeloid
  marker_genes:
    Adipocyte: [ADIPOQ, GPAM]
    Myeloid: [MRC1, LPL]
myeloid:
  fine_celltype_key: cell_state_am
  myeloid_labels: [Myeloid, Macrophage]
  condition_key: condition
  color_key: cell_state_am
compass:
  cell_state_h5ad: data/swat_all_seurat_integration_scott_ref_stringent_rna_global_annotated_cell_states.h5ad
  sample_metadata: metadata/macro_samples.tsv
  flux_matrix: results/compass/macrophage_fluxes.tsv
  expression_output: results/compass/expression_matrix_log1p_per_sample_mean.tsv
  condition_key: condition
  condition_reference: LN
  condition_test: OB
  condition_alternative: WL
lam_trm:
  flux_matrix: results/compass/macrophage_fluxes.tsv
  sample_metadata: metadata/macro_samples.tsv
  subtype_key: macrophage_subtype
  lam_labels: [LAM]
  trm_labels: [TRM]
scenith:
  scenith_csv: data/scenith_ob_macrophages.csv
```

Each section below documents the analysis logic encapsulated by the runner for transparency and reproducibility.

---

## 1. Shared pre-processing and broad clustering

Source script: [`single-nuc/clustering/global_clustering_and_integration.py`](../clustering/global_clustering_and_integration.py)

1. **Load the integrated dataset**
   - Read the integrated AnnData object and confirm the presence of raw counts in `adata.layers["raw_counts"]`.
   - Apply `sc.pp.normalize_total` (target sum 10,000) followed by `sc.pp.log1p`, storing results in `adata.layers["log1p_counts"]`.

2. **Regress technical covariates and compute latent space**
   - Regress mitochondrial content (`mt.percent`), ribosomal content (`ribo.percent`), and total counts (`nCount_RNA`) using `sc.pp.regress_out`.
   - Identify highly variable genes (`sc.pp.highly_variable_genes`) and compute a 40-component PCA (`sc.tl.pca`).
   - Apply Harmony batch correction with `sample` as the batch key, then build a bbknn neighbour graph on the Harmony embedding.
   - Generate UMAP coordinates (`sc.tl.umap`) and run Leiden clustering at low resolution to obtain major cell populations.

3. **Annotate broad cell types**
   - Use marker genes defined in the script (e.g., *ADIPOQ*, *GPAM*, *VWF*, *CDH5*, *MRC1*) to create dot plots and assign each Leiden cluster to a broad cell type.
   - Replace numeric Leiden labels with human-readable cell type names and save the annotated object with the field `adata.obs.cell_type_am`.

4. **Save per-cell-type AnnData files**
   - Write the annotated global object to `*_global_annotated.h5ad`.
   - Iterate over `adata.obs.cell_type_am` and export one AnnData file per broad type (e.g., `*_global_annotated_Myeloid.h5ad`). These serve as inputs for cell-type-specific analyses below.

---

## 2. Figure 2a – Myeloid UMAP and density plots

Inputs: `*_global_annotated_Myeloid.h5ad`

1. **Subset myeloid classes**
   - Restrict the object to myeloid subclasses using `cell_type_am_fine` or `cell_state_am` annotations.

2. **Compute myeloid UMAP**
   - Ensure counts are normalised/log-transformed (`adata.layers["log1p_counts"]`).
   - Run `sc.pp.neighbors` and `sc.tl.umap` on the myeloid subset.
   - Visualise the embedding coloured by myeloid class using `sc.pl.umap` (top panel of Fig. 2a).

3. **Estimate density per condition**
   - Apply `sc.tl.embedding_density` using `condition` (LN, OB, WL) to compute density values in UMAP space.
   - Plot condition-specific density maps with `sc.pl.embedding_density` (bottom panels of Fig. 2a).

4. **Export results**
   - Save UMAP coordinates and density estimates to TSV/CSV files for reproducibility.
   - Store high-resolution images for the class-coloured UMAP and density panels.

---

## 3. Figure 2e – Compass-based metabolic activation in macrophages

Inputs: `*_global_annotated_cell_states.h5ad`

1. **Export per-sample expression matrices**
   - Use [`single-nuc/compass/compass_matrix_export.py`](../compass/compass_matrix_export.py) with `pool_key = "sample"`.
   - Set `adata.X = adata.layers["log1p_counts"]` and subset to macrophage cells (`cell_type_am_fine == "Macrophage"`).
   - Call `get_df` to compute mean log1p expression per gene per sample; export as a tab-delimited file for Compass.

2. **Run Compass**
   - Execute [`single-nuc/compass/run_compass.sh`](../compass/run_compass.sh) with the exported matrix via the Compass CLI (`--species homo_sapiens`).
   - Capture reaction-level flux estimates for all samples.

3. **Compute condition effect sizes**
   - Load Compass outputs and group samples by `condition` (LN, OB, WL).
   - For each reaction, calculate Cohen's *d* for OB vs LN and OB vs WL, run Wilcoxon rank-sum tests, and apply FDR correction.
   - Label reactions as OB-high (positive *d*, FDR < 0.05), OB-low (negative *d*, FDR < 0.05), or non-significant, and summarise counts per metabolic pathway (e.g., PPP, OXPHOS, glycolysis/gluconeogenesis, fatty-acid pathways).

4. **Visualise global and pathway activation**
   - Produce a heat/strip plot of reaction effect sizes coloured by activation class.
   - Create pathway-level summary plots (bar or pie charts) that report the proportion of significant reactions, noting sample sizes (n = 24 LN, 25 OB, 24 WL).

---

## 4. Figure 2g – LAM vs TRM metabolic comparison and SCENITH

1. **LAM vs TRM Compass analysis**
   - Annotate macrophage Compass samples with subtype labels (LAM vs TRM) using cell-state metadata.
   - Compute mean flux per subtype, calculate Cohen's *d*, run FDR correction, and classify reactions as LAM-high, LAM-low, or non-significant (sample sizes: MYE1 n = 86, MYE2 n = 74, MYE3 n = 80).
   - Visualise subtype differences with a heat/strip plot similar to Figure 2e.

2. **SCENITH assay analysis**
   - Import SCENITH flow cytometry data measuring basal respiration (HPG incorporation) and glycolytic capacity (delta HPG) for LAM and TRM populations from OB donors (n = 7).
   - Compute donor-level means, perform paired *t*-tests, and plot paired dot or box plots showing LAM vs TRM values with mean ± s.e.m.
   - These plots form the lower panel of Figure 2g.

---

## Outputs summary

| Figure panel | Key deliverables |
|--------------|------------------|
| Fig. 2a      | Myeloid UMAP embedding, per-condition density plots, exported coordinates/density tables |
| Fig. 2e      | Compass flux matrices, effect-size tables, global and pathway activation plots |
| Fig. 2g      | LAM vs TRM Compass comparison plots, SCENITH paired analyses |

This document consolidates the existing scripts and analysis steps required to reproduce the single-nucleus components of Figure 2.
