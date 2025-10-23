from __future__ import annotations

import json
import math
import subprocess
from dataclasses import dataclass, field, fields
from pathlib import Path
from collections.abc import Mapping as MappingABC, Sequence as SequenceABC
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union, get_args, get_origin, TYPE_CHECKING

try:  # third-party scientific stack (optional in mock mode)
    import matplotlib.pyplot as _plt  # type: ignore
except ImportError:  # pragma: no cover - optional dependency
    _plt = None

try:  # pragma: no cover - optional dependency
    import numpy as _np  # type: ignore
except ImportError:
    _np = None

try:  # pragma: no cover - optional dependency
    import pandas as _pd  # type: ignore
except ImportError:
    _pd = None

try:  # pragma: no cover - optional dependency
    import scanpy as _sc  # type: ignore
except ImportError:
    _sc = None

try:  # pragma: no cover - optional dependency
    from scipy import stats as _stats  # type: ignore
except ImportError:
    _stats = None

try:  # optional dependency but widely available
    import yaml
except ImportError:  # pragma: no cover - fall back to json only
    yaml = None

if TYPE_CHECKING:  # pragma: no cover - type checking only
    from anndata import AnnData
    import pandas as pd
    import numpy as np


def _require_dependency(name: str, module):
    if module is None:
        raise RuntimeError(
            f"Optional dependency '{name}' is required for the full Figure 2 pipeline. "
            "Install the scientific stack or set 'mock_mode: true' in the configuration to "
            "generate placeholder outputs without third-party packages."
        )
    return module


@dataclass
class PreprocessingConfig:
    """Configuration for the global preprocessing workflow."""

    regress_covariates: Sequence[str] = ("mt.percent", "ribo.percent", "nCount_RNA")
    n_top_genes: int = 4000
    n_pcs: int = 40
    harmony_batch_key: str = "sample"
    bbknn_neighbors: int = 10
    leiden_resolution: float = 0.25
    marker_genes: Mapping[str, Sequence[str]] = field(default_factory=dict)
    cluster_cell_type_map: Mapping[str, str] = field(default_factory=dict)
    cell_type_column: str = "cell_type_am"
    raw_layer: str = "raw_counts"
    log1p_layer: str = "log1p_counts"
    global_output_name: str = "global_annotated.h5ad"
    split_cell_type_column: Optional[str] = "cell_type_am"


@dataclass
class MyeloidConfig:
    """Settings for reproducing the myeloid UMAP and density plots."""

    enabled: bool = True
    fine_celltype_key: str = "cell_type_am_fine"
    myeloid_labels: Sequence[str] = field(default_factory=lambda: ["Myeloid", "Macrophage", "Monocyte"])
    umap_neighbors: int = 15
    umap_min_dist: float = 0.3
    umap_spread: float = 1.0
    condition_key: str = "condition"
    color_key: str = "cell_state_am"
    output_prefix: str = "figure2a"


@dataclass
class CompassConfig:
    """Configuration for Compass flux analyses (Figure 2e)."""

    enabled: bool = True
    cell_state_h5ad: Optional[Path] = None
    sample_key: str = "sample"
    layer: str = "log1p_counts"
    cell_type_key: str = "cell_type_am_fine"
    macrophage_labels: Sequence[str] = field(default_factory=lambda: ["Macrophage"])
    expression_output: Path = Path("expression_matrix_log1p_per_sample_mean.tsv")
    compass_command: Optional[Sequence[str]] = None
    compass_workdir: Optional[Path] = None
    flux_matrix: Optional[Path] = None
    sample_metadata: Optional[Path] = None
    condition_key: str = "condition"
    condition_reference: str = "LN"
    condition_test: str = "OB"
    condition_alternative: Optional[str] = "WL"
    pathway_annotation: Optional[Path] = None
    pathway_column: str = "pathway"
    reaction_column: str = "reaction"
    significance_alpha: float = 0.05
    output_prefix: str = "figure2e"

    def __post_init__(self) -> None:
        if self.cell_state_h5ad is not None:
            self.cell_state_h5ad = Path(self.cell_state_h5ad)
        if self.expression_output:
            self.expression_output = Path(self.expression_output)
        if self.compass_workdir is not None:
            self.compass_workdir = Path(self.compass_workdir)
        if self.flux_matrix is not None:
            self.flux_matrix = Path(self.flux_matrix)
        if self.sample_metadata is not None:
            self.sample_metadata = Path(self.sample_metadata)
        if self.pathway_annotation is not None:
            self.pathway_annotation = Path(self.pathway_annotation)


@dataclass
class LamTrmConfig:
    """Configuration for LAM vs TRM comparison (Figure 2g upper panel)."""

    enabled: bool = True
    flux_matrix: Optional[Path] = None
    sample_metadata: Optional[Path] = None
    subtype_key: str = "macrophage_subtype"
    lam_labels: Sequence[str] = field(default_factory=lambda: ["LAM"])
    trm_labels: Sequence[str] = field(default_factory=lambda: ["TRM"])
    pathway_annotation: Optional[Path] = None
    pathway_column: str = "pathway"
    reaction_column: str = "reaction"
    significance_alpha: float = 0.05
    output_prefix: str = "figure2g_flux"

    def __post_init__(self) -> None:
        if self.flux_matrix is not None:
            self.flux_matrix = Path(self.flux_matrix)
        if self.sample_metadata is not None:
            self.sample_metadata = Path(self.sample_metadata)
        if self.pathway_annotation is not None:
            self.pathway_annotation = Path(self.pathway_annotation)


@dataclass
class ScenithConfig:
    """Configuration for the SCENITH paired analysis (Figure 2g lower panel)."""

    enabled: bool = True
    scenith_csv: Optional[Path] = None
    donor_column: str = "donor"
    subtype_column: str = "subtype"
    lam_label: str = "LAM"
    trm_label: str = "TRM"
    basal_column: str = "basal_respiration"
    glycolysis_column: str = "glycolytic_capacity"
    output_prefix: str = "figure2g_scenith"

    def __post_init__(self) -> None:
        if self.scenith_csv is not None:
            self.scenith_csv = Path(self.scenith_csv)


@dataclass
class Figure2Config:
    """Top-level configuration for the Figure 2 pipeline."""

    integrated_h5ad: Path
    output_dir: Path
    preprocess: PreprocessingConfig = field(default_factory=PreprocessingConfig)
    myeloid: MyeloidConfig = field(default_factory=MyeloidConfig)
    compass: CompassConfig = field(default_factory=CompassConfig)
    lam_trm: LamTrmConfig = field(default_factory=LamTrmConfig)
    scenith: ScenithConfig = field(default_factory=ScenithConfig)
    mock_mode: bool = False

    @classmethod
    def from_file(cls, path: Path) -> "Figure2Config":
        data = _read_structured_file(path)
        return _build_dataclass(cls, data)

    def prepare_output_dirs(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        (self.output_dir / "logs").mkdir(exist_ok=True)

    def __post_init__(self) -> None:
        self.integrated_h5ad = Path(self.integrated_h5ad)
        self.output_dir = Path(self.output_dir)


def _read_structured_file(path: Path) -> Mapping[str, object]:
    with Path(path).open("r", encoding="utf-8") as handle:
        text = handle.read()
    if path.suffix.lower() in {".json"}:
        return json.loads(text)
    if path.suffix.lower() in {".yml", ".yaml"}:
        if yaml is None:
            raise RuntimeError("PyYAML is required to parse YAML configuration files")
        return yaml.safe_load(text)
    raise ValueError(f"Unsupported config extension for {path}")


def _build_dataclass(cls, data: Mapping[str, object]):
    if not isinstance(data, Mapping):
        raise TypeError(f"Expected mapping to build {cls.__name__}")
    init_kwargs = {}
    for field_def in fields(cls):
        name = field_def.name
        if name not in data:
            continue
        value = data[name]
        init_kwargs[name] = _coerce_value(field_def.type, value)
    return cls(**init_kwargs)  # type: ignore[arg-type]


def _coerce_value(field_type, value):
    if isinstance(field_type, str):
        field_type = globals().get(field_type, field_type)
    origin = get_origin(field_type)
    if origin is None:
        if hasattr(field_type, "__dataclass_fields__"):
            return _build_dataclass(field_type, value)
        if field_type is Path and not isinstance(value, Path):
            return Path(value)
        return value
    if origin in {list, tuple, set} or (origin is not None and isinstance(origin, type) and issubclass(origin, SequenceABC)):
        args = get_args(field_type)
        if not isinstance(value, (list, tuple, set)):
            return value
        container_type = list if origin is list else type(value)
        if args:
            return container_type(_coerce_value(args[0], item) for item in value)
        return container_type(value)
    if origin is dict or (origin is not None and isinstance(origin, type) and issubclass(origin, MappingABC)):
        args = get_args(field_type)
        key_type = args[0] if len(args) > 0 else None
        val_type = args[1] if len(args) > 1 else None
        return {
            _coerce_value(key_type, key) if key_type else key:
            _coerce_value(val_type, val) if val_type else val
            for key, val in value.items()
        }
    if origin is Union:
        for arg in get_args(field_type):
            if arg is type(None):  # noqa: E721
                if value is None:
                    return None
                continue
            try:
                return _coerce_value(arg, value)
            except Exception:
                continue
        return value
    return value


# ---------------------------------------------------------------------------
# Pipeline orchestration
# ---------------------------------------------------------------------------

def run_pipeline(config: Figure2Config) -> None:
    config.prepare_output_dirs()

    if config.mock_mode:
        run_mock_pipeline(config)
        return

    sc = _require_dependency("scanpy", _sc)
    adata = sc.read_h5ad(config.integrated_h5ad)
    print(f"Loaded integrated object with {adata.n_obs} cells and {adata.n_vars} genes")

    global_adata = run_global_preprocessing(adata, config)
    global_path = config.output_dir / config.preprocess.global_output_name
    print(f"Saving globally annotated object to {global_path}")
    global_adata.write_h5ad(global_path)

    if config.preprocess.split_cell_type_column:
        write_cell_type_splits(global_adata, config)

    if config.myeloid.enabled:
        run_myeloid_umap(global_adata, config)

    if config.compass.enabled:
        run_compass_workflow(config)

    if config.lam_trm.enabled:
        run_lam_trm_workflow(config)

    if config.scenith.enabled:
        run_scenith_workflow(config)


def run_mock_pipeline(config: Figure2Config) -> None:
    print("Running Figure 2 pipeline in mock mode; generating placeholder artefacts.")
    base = config.output_dir
    log_path = base / "logs" / "mock_run.log"
    log_path.write_text(
        "Mock mode executed. No scientific dependencies were required; files are placeholders.\n",
        encoding="utf-8",
    )

    _write_placeholder(
        base / config.preprocess.global_output_name,
        "Mock AnnData content: placeholder for global annotated object.\n",
    )

    # Myeloid outputs
    myeloid_prefix = config.myeloid.output_prefix
    _write_pdf_placeholder(base / f"{myeloid_prefix}_classes.pdf", "Mock myeloid UMAP")
    _write_pdf_placeholder(base / f"{myeloid_prefix}_density_mock.pdf", "Mock density plot")
    coords_path = base / f"{myeloid_prefix}_umap_coordinates.tsv"
    _write_tsv(
        coords_path,
        ["cell_id", "UMAP1", "UMAP2", config.myeloid.condition_key, config.myeloid.color_key],
        [
            ["cell_1", "0.1", "0.2", "LN", "MockMyeloid"],
            ["cell_2", "0.2", "0.3", "OB", "MockMyeloid"],
            ["cell_3", "0.3", "0.1", "WL", "MockMyeloid"],
        ],
    )

    # Compass placeholder artefacts
    compass = config.compass
    expression_path = _resolve_output_path(base, compass.expression_output)
    expression_path.parent.mkdir(parents=True, exist_ok=True)
    _write_tsv(
        expression_path,
        ["gene", "Sample_LN", "Sample_OB", "Sample_WL"],
        [["G1", "1.0", "2.0", "1.5"], ["G2", "0.5", "1.2", "0.8"]],
    )

    metadata_path = None
    if compass.sample_metadata is not None:
        metadata_path = _resolve_output_path(base, compass.sample_metadata)
        metadata_path.parent.mkdir(parents=True, exist_ok=True)
        _write_tsv(
            metadata_path,
            [compass.sample_key, compass.condition_key, config.lam_trm.subtype_key],
            [
                ["Sample_LN", compass.condition_reference, "TRM"],
                ["Sample_OB", compass.condition_test, "LAM"],
                ["Sample_WL", compass.condition_alternative or "WL", "LAM"],
            ],
        )

    flux_path = None
    if compass.flux_matrix is not None:
        flux_path = _resolve_output_path(base, compass.flux_matrix)
    else:
        flux_path = base / "mock_flux_matrix.tsv"
    flux_path.parent.mkdir(parents=True, exist_ok=True)
    _write_tsv(
        flux_path,
        ["reaction", "Sample_LN", "Sample_OB", "Sample_WL"],
        [["R1", "0.5", "1.0", "0.7"], ["R2", "0.2", "0.6", "0.4"], ["R3", "0.8", "0.4", "0.5"]],
    )

    _write_pdf_placeholder(
        base / f"{compass.output_prefix}_{compass.condition_test}_vs_{compass.condition_reference}.pdf",
        "Mock Compass comparison",
    )
    if compass.condition_alternative:
        _write_pdf_placeholder(
            base / f"{compass.output_prefix}_{compass.condition_test}_vs_{compass.condition_alternative}.pdf",
            "Mock Compass alternative comparison",
        )

    # LAM/TRM artefact (re-uses flux metadata)
    lam_cfg = config.lam_trm
    if lam_cfg.sample_metadata and metadata_path is None:
        metadata_path = _resolve_output_path(base, lam_cfg.sample_metadata)
        metadata_path.parent.mkdir(parents=True, exist_ok=True)
        _write_tsv(
            metadata_path,
            [config.compass.sample_key, lam_cfg.subtype_key],
            [["Sample_OB", "LAM"], ["Sample_LN", "TRM"]],
        )
    _write_pdf_placeholder(base / f"{lam_cfg.output_prefix}.pdf", "Mock LAM vs TRM plot")

    # SCENITH placeholders
    scenith_cfg = config.scenith
    scenith_path = _resolve_output_path(base, scenith_cfg.scenith_csv) if scenith_cfg.scenith_csv else base / "mock_scenith.csv"
    scenith_path.parent.mkdir(parents=True, exist_ok=True)
    _write_tsv(
        scenith_path,
        [
            scenith_cfg.donor_column,
            scenith_cfg.subtype_column,
            scenith_cfg.basal_column,
            scenith_cfg.glycolysis_column,
        ],
        [
            ["Donor1", scenith_cfg.lam_label, "1.2", "0.8"],
            ["Donor1", scenith_cfg.trm_label, "0.9", "0.6"],
        ],
    )
    _write_pdf_placeholder(base / f"{scenith_cfg.output_prefix}.pdf", "Mock SCENITH plot")

    print(f"Mock artefacts written to {base}")

# ---------------------------------------------------------------------------
# Global preprocessing
# ---------------------------------------------------------------------------

def run_global_preprocessing(adata: 'AnnData', config: Figure2Config) -> 'AnnData':
    sc = _require_dependency("scanpy", _sc)
    cfg = config.preprocess
    if cfg.raw_layer not in adata.layers:
        adata.layers[cfg.raw_layer] = adata.X.copy()
    if cfg.log1p_layer not in adata.layers:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.layers[cfg.log1p_layer] = adata.X.copy()
    else:
        adata.X = adata.layers[cfg.log1p_layer].copy()

    regressors = list(cfg.regress_covariates)
    available_regressors = [cov for cov in regressors if cov in adata.obs]
    if available_regressors:
        print(f"Regressing out covariates: {available_regressors}")
        sc.pp.regress_out(adata, keys=available_regressors)
    else:
        print("No regression covariates present; skipping regress_out")

    sc.pp.highly_variable_genes(adata, n_top_genes=cfg.n_top_genes)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=cfg.n_pcs)

    sc.external.pp.harmony_integrate(adata, key=cfg.harmony_batch_key, basis="X_pca")
    sc.external.pp.bbknn(
        adata,
        batch_key=cfg.harmony_batch_key,
        neighbors_within_batch=cfg.bbknn_neighbors,
        use_rep="X_pca_harmony",
    )
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=cfg.leiden_resolution)

    annotate_clusters(adata, config)
    if cfg.marker_genes:
        save_marker_dotplot(adata, cfg, config.output_dir)

    return adata


def annotate_clusters(adata: 'AnnData', config: Figure2Config) -> None:
    cfg = config.preprocess
    if not cfg.cluster_cell_type_map:
        print("No cluster-to-cell-type mapping provided; retaining Leiden labels")
        adata.obs[cfg.cell_type_column] = adata.obs["leiden"].astype(str)
        return

    cell_type_labels = adata.obs["leiden"].astype(str).map(cfg.cluster_cell_type_map)
    adata.obs[cfg.cell_type_column] = cell_type_labels.fillna(adata.obs["leiden"].astype(str))


def save_marker_dotplot(adata: 'AnnData', cfg: PreprocessingConfig, output_dir: Path) -> None:
    sc = _require_dependency("scanpy", _sc)
    plt = _require_dependency("matplotlib", _plt)
    marker_genes: Dict[str, Sequence[str]] = dict(cfg.marker_genes)
    genes = sorted({gene for gene_list in marker_genes.values() for gene in gene_list})
    genes_present = [gene for gene in genes if gene in adata.var_names]
    if not genes_present:
        print("None of the requested marker genes were found; skipping dotplot")
        return
    categories = list(marker_genes.keys())
    marker_matrix = {cat: [gene for gene in marker_genes[cat] if gene in genes_present] for cat in categories}
    plot_obj = sc.pl.dotplot(
        adata,
        var_names=marker_matrix,
        groupby="leiden",
        standard_scale="var",
        show=False,
    )
    dotplot_path = output_dir / "figure2_global_marker_dotplot.pdf"
    if hasattr(plot_obj, "savefig"):
        plot_obj.savefig(dotplot_path)
        if hasattr(plot_obj, "figure") and plot_obj.figure is not None:
            plt.close(plot_obj.figure)
    else:
        fig = plot_obj if isinstance(plot_obj, plt.Figure) else plot_obj.figure
        fig.savefig(dotplot_path, bbox_inches="tight")
        plt.close(fig)
    print(f"Saved marker gene dotplot to {dotplot_path}")


def write_cell_type_splits(adata: 'AnnData', config: Figure2Config) -> None:
    column = config.preprocess.split_cell_type_column
    if column is None or column not in adata.obs:
        print("Requested split column missing; skipping per-cell-type exports")
        return
    base_name = Path(config.preprocess.global_output_name).with_suffix("")
    for cell_type, subset in adata.obs.groupby(column):
        if not isinstance(cell_type, str):
            cell_type = str(cell_type)
        safe_name = cell_type.replace("/", "-").replace(" ", "_")
        path = config.output_dir / f"{base_name}_{safe_name}.h5ad"
        adata[subset.index].write_h5ad(path)
        print(f"Wrote subset for {cell_type} with {subset.index.size} cells -> {path}")


# ---------------------------------------------------------------------------
# Figure 2a - Myeloid UMAP and density
# ---------------------------------------------------------------------------

def run_myeloid_umap(global_adata: 'AnnData', config: Figure2Config) -> None:
    sc = _require_dependency("scanpy", _sc)
    plt = _require_dependency("matplotlib", _plt)
    pd = _require_dependency("pandas", _pd)
    cfg = config.myeloid
    color_key = cfg.color_key if cfg.color_key in global_adata.obs.columns else "leiden"
    if cfg.fine_celltype_key not in global_adata.obs:
        raise KeyError(f"Column '{cfg.fine_celltype_key}' not found in AnnData.obs")
    myeloid_mask = global_adata.obs[cfg.fine_celltype_key].isin(cfg.myeloid_labels)
    if myeloid_mask.sum() == 0:
        raise ValueError("No myeloid cells found with the specified labels")
    adata = global_adata[myeloid_mask].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000)
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=30)
    sc.pp.neighbors(adata, n_neighbors=cfg.umap_neighbors, use_rep="X_pca")
    sc.tl.umap(adata, min_dist=cfg.umap_min_dist, spread=cfg.umap_spread)

    plot = sc.pl.umap(
        adata,
        color=color_key,
        legend_loc="on data",
        title="Myeloid subclasses",
        show=False,
    )
    fig = plot if isinstance(plot, plt.Figure) else plot.figure
    umap_path = config.output_dir / f"{cfg.output_prefix}_classes.pdf"
    fig.savefig(umap_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved myeloid class UMAP to {umap_path}")

    if cfg.condition_key not in adata.obs:
        raise ValueError(f"Condition key '{cfg.condition_key}' not present for density plot")
    sc.tl.embedding_density(adata, basis="umap", groupby=cfg.condition_key)
    conditions = adata.obs[cfg.condition_key].unique()
    for condition in conditions:
        density_plot = sc.pl.embedding_density(
            adata,
            basis="umap",
            key=condition,
            title=f"{condition} density",
            show=False,
        )
        density_fig = density_plot if isinstance(density_plot, plt.Figure) else density_plot.figure
        density_path = config.output_dir / f"{cfg.output_prefix}_density_{condition}.pdf"
        density_fig.savefig(density_path, bbox_inches="tight")
        plt.close(density_fig)
        print(f"Saved myeloid density map for {condition} to {density_path}")

    coords = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"], index=adata.obs_names)
    coords[cfg.condition_key] = adata.obs[cfg.condition_key].values
    coords[color_key] = adata.obs[color_key].values
    coords_path = config.output_dir / f"{cfg.output_prefix}_umap_coordinates.tsv"
    coords.to_csv(coords_path, sep="\t")
    print(f"Saved UMAP coordinates to {coords_path}")


# ---------------------------------------------------------------------------
# Figure 2e - Compass flux analysis
# ---------------------------------------------------------------------------

def run_compass_workflow(config: Figure2Config) -> None:
    cfg = config.compass
    if cfg.cell_state_h5ad is None:
        raise ValueError("Compass configuration requires 'cell_state_h5ad'")

    sc = _require_dependency("scanpy", _sc)
    pd = _require_dependency("pandas", _pd)
    plt = _require_dependency("matplotlib", _plt)
    adata = sc.read_h5ad(cfg.cell_state_h5ad)
    print(f"Loaded cell-state object with {adata.n_obs} cells for Compass export")

    if cfg.layer:
        if cfg.layer not in adata.layers:
            raise KeyError(f"Layer '{cfg.layer}' not found in AnnData layers")
        adata.X = adata.layers[cfg.layer].copy()
    if adata.raw is None:
        adata.raw = adata

    macrophage_mask = adata.obs[cfg.cell_type_key].isin(cfg.macrophage_labels)
    subset = adata[macrophage_mask].copy()
    if subset.n_obs == 0:
        raise ValueError("No macrophage cells found for Compass export")

    expression_matrix = compute_mean_expression(subset, cfg.sample_key)
    expression_path = Path(cfg.expression_output)
    if not expression_path.is_absolute():
        expression_path = config.output_dir / expression_path
    expression_path.parent.mkdir(parents=True, exist_ok=True)
    expression_matrix.to_csv(expression_path, sep="\t")
    print(f"Exported Compass input matrix to {expression_path}")

    if cfg.compass_command:
        workdir = Path(cfg.compass_workdir or expression_path.parent)
        run_subprocess(cfg.compass_command, cwd=workdir)

    if cfg.flux_matrix:
        flux_path = Path(cfg.flux_matrix)
    else:
        flux_path = expression_path.parent / "fluxes.tsv"
    if not flux_path.exists():
        print(f"Flux matrix {flux_path} not found; skipping downstream Compass plots")
        return

    if cfg.sample_metadata is None:
        raise ValueError("Compass downstream analysis requires 'sample_metadata'")
    metadata = pd.read_csv(cfg.sample_metadata, sep=None, engine="python")
    flux = load_flux_matrix(flux_path)

    merged = metadata.set_index(cfg.sample_key).join(flux.T, how="inner")
    group_columns = [cfg.condition_reference, cfg.condition_test]
    if cfg.condition_alternative:
        group_columns.append(cfg.condition_alternative)
    available_groups = merged[cfg.condition_key].unique()
    print(f"Available conditions for Compass analysis: {available_groups}")

    fig = plot_condition_effects(
        flux,
        merged,
        group_a=cfg.condition_reference,
        group_b=cfg.condition_test,
        condition_key=cfg.condition_key,
        pathway_annotation=cfg.pathway_annotation,
        pathway_column=cfg.pathway_column,
        reaction_column=cfg.reaction_column,
        alpha=cfg.significance_alpha,
        title=f"{cfg.condition_test} vs {cfg.condition_reference}"
    )
    output_path = config.output_dir / f"{cfg.output_prefix}_{cfg.condition_test}_vs_{cfg.condition_reference}.pdf"
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved Compass comparison plot to {output_path}")

    if cfg.condition_alternative:
        fig = plot_condition_effects(
            flux,
            merged,
            group_a=cfg.condition_alternative,
            group_b=cfg.condition_test,
            condition_key=cfg.condition_key,
            pathway_annotation=cfg.pathway_annotation,
            pathway_column=cfg.pathway_column,
            reaction_column=cfg.reaction_column,
            alpha=cfg.significance_alpha,
            title=f"{cfg.condition_test} vs {cfg.condition_alternative}"
        )
        output_path = config.output_dir / f"{cfg.output_prefix}_{cfg.condition_test}_vs_{cfg.condition_alternative}.pdf"
        fig.savefig(output_path, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved Compass comparison plot to {output_path}")


# ---------------------------------------------------------------------------
# Figure 2g - LAM vs TRM flux comparison
# ---------------------------------------------------------------------------

def run_lam_trm_workflow(config: Figure2Config) -> None:
    cfg = config.lam_trm
    if not cfg.flux_matrix or not Path(cfg.flux_matrix).exists():
        raise ValueError("LAM/TRM configuration requires an existing flux matrix")
    if cfg.sample_metadata is None:
        raise ValueError("LAM/TRM configuration requires sample metadata")

    pd = _require_dependency("pandas", _pd)
    plt = _require_dependency("matplotlib", _plt)
    flux = load_flux_matrix(cfg.flux_matrix)
    metadata = pd.read_csv(cfg.sample_metadata, sep=None, engine="python").set_index(config.compass.sample_key)
    merged = metadata.join(flux.T, how="inner")

    fig = plot_condition_effects(
        flux,
        merged,
        group_a=cfg.lam_labels,
        group_b=cfg.trm_labels,
        condition_key=cfg.subtype_key,
        pathway_annotation=cfg.pathway_annotation,
        pathway_column=cfg.pathway_column,
        reaction_column=cfg.reaction_column,
        alpha=cfg.significance_alpha,
        title="LAM vs TRM",
    )
    output_path = config.output_dir / f"{cfg.output_prefix}.pdf"
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved LAM vs TRM flux plot to {output_path}")


# ---------------------------------------------------------------------------
# Figure 2g - SCENITH paired analysis
# ---------------------------------------------------------------------------

def run_scenith_workflow(config: Figure2Config) -> None:
    cfg = config.scenith
    if not cfg.scenith_csv or not Path(cfg.scenith_csv).exists():
        raise ValueError("SCENITH configuration requires an input CSV file")

    pd = _require_dependency("pandas", _pd)
    plt = _require_dependency("matplotlib", _plt)
    data = pd.read_csv(cfg.scenith_csv, sep=None, engine="python")
    if cfg.subtype_column not in data.columns:
        raise KeyError(f"Missing subtype column '{cfg.subtype_column}' in SCENITH data")

    subset = data[data[cfg.subtype_column].isin([cfg.lam_label, cfg.trm_label])]
    stats_summary = compute_paired_statistics(
        subset,
        donor_key=cfg.donor_column,
        subtype_key=cfg.subtype_column,
        lam_label=cfg.lam_label,
        trm_label=cfg.trm_label,
        value_columns=[cfg.basal_column, cfg.glycolysis_column],
    )
    plot_scenith_pairs(
        subset,
        cfg,
        stats_summary,
        output_dir=config.output_dir,
    )


# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def compute_mean_expression(adata: 'AnnData', groupby: str) -> 'pd.DataFrame':
    np = _require_dependency("numpy", _np)
    pd = _require_dependency("pandas", _pd)
    if groupby not in adata.obs:
        raise KeyError(f"Grouping column '{groupby}' is missing")
    matrices = []
    sample_names = []
    for sample, df in adata.obs.groupby(groupby):
        sample_names.append(sample)
        matrices.append(np.asarray(adata[df.index].X.mean(axis=0)).ravel())
    matrix = np.vstack(matrices).T
    return pd.DataFrame(matrix, index=adata.var_names, columns=sample_names)


def load_flux_matrix(path: Path | str) -> 'pd.DataFrame':
    pd = _require_dependency("pandas", _pd)
    df = pd.read_csv(path, sep=None, engine="python")
    if df.shape[1] <= df.shape[0]:
        df = df.set_index(df.columns[0])
    return df


def run_subprocess(command: Sequence[str], cwd: Optional[Path] = None) -> None:
    print(f"Running command: {' '.join(command)} (cwd={cwd})")
    subprocess.run(command, cwd=cwd, check=True)


def plot_condition_effects(
    flux: 'pd.DataFrame',
    metadata: 'pd.DataFrame',
    group_a: Sequence[str] | str,
    group_b: Sequence[str] | str,
    condition_key: str,
    pathway_annotation: Optional[Path],
    pathway_column: str,
    reaction_column: str,
    alpha: float,
    title: str,
) -> object:
    np = _require_dependency("numpy", _np)
    pd = _require_dependency("pandas", _pd)
    plt = _require_dependency("matplotlib", _plt)
    stats = _require_dependency("scipy.stats", _stats)

    if isinstance(group_a, str):
        group_a = [group_a]
    if isinstance(group_b, str):
        group_b = [group_b]
    mask_a = metadata[condition_key].isin(group_a)
    mask_b = metadata[condition_key].isin(group_b)
    group_a_samples = [sample for sample in metadata.index[mask_a] if sample in flux.columns]
    group_b_samples = [sample for sample in metadata.index[mask_b] if sample in flux.columns]
    if len(group_a_samples) == 0 or len(group_b_samples) == 0:
        raise ValueError("One of the comparison groups has no samples")

    effect_records: List[Dict[str, object]] = []
    for reaction in flux.index:
        values_a = flux.loc[reaction, group_a_samples].dropna().to_numpy(dtype=float)
        values_b = flux.loc[reaction, group_b_samples].dropna().to_numpy(dtype=float)
        if values_a.size < 2 or values_b.size < 2:
            continue
        effect = cohen_d(values_b, values_a)
        stat, p_value = stats.ranksums(values_b, values_a)
        record = {
            "reaction": reaction,
            "effect_size": effect,
            "p_value": p_value,
        }
        effect_records.append(record)

    results = pd.DataFrame(effect_records)
    if results.empty:
        raise ValueError("No reactions available after filtering for effect sizes")
    results["fdr"] = _benjamini_hochberg(results["p_value"].values)
    if pathway_annotation:
        pathway_df = pd.read_csv(pathway_annotation, sep=None, engine="python")
        results = results.merge(
            pathway_df[[reaction_column, pathway_column]],
            left_on="reaction",
            right_on=reaction_column,
            how="left",
        )
    else:
        results[pathway_column] = "Unassigned"

    results["classification"] = np.where(
        (results["effect_size"] > 0) & (results["fdr"] < alpha),
        "higher",
        np.where((results["effect_size"] < 0) & (results["fdr"] < alpha), "lower", "ns"),
    )

    pivot = (
        results.groupby(pathway_column)["classification"]
        .value_counts(normalize=True)
        .rename("fraction")
        .reset_index()
    )
    fig, axes = plt.subplots(2, 1, figsize=(8, 10), gridspec_kw={"height_ratios": [3, 1]})
    ax = axes[0]
    ax.scatter(
        results["effect_size"],
        -np.log10(results["fdr"].replace(0, np.nextafter(0, 1))),
        c=results["classification"].map({"higher": "#b2182b", "lower": "#2166ac", "ns": "#cccccc"}),
        alpha=0.7,
    )
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("Effect size (Cohen's d)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title(title)

    ax2 = axes[1]
    fraction_table = pivot.pivot(index=pathway_column, columns="classification", values="fraction").fillna(0)
    fraction_table.loc[:, [col for col in ["higher", "lower", "ns"] if col in fraction_table.columns]].plot(
        kind="bar",
        stacked=True,
        color={"higher": "#b2182b", "lower": "#2166ac", "ns": "#cccccc"},
        ax=ax2,
    )
    ax2.set_ylabel("Fraction of reactions")
    ax2.set_xlabel("Pathway")
    ax2.legend(title="Classification", bbox_to_anchor=(1.05, 1), loc="upper left")
    fig.tight_layout()
    return fig


def cohen_d(group1, group2) -> float:
    np = _require_dependency("numpy", _np)
    diff = group1.mean() - group2.mean()
    n1 = group1.size
    n2 = group2.size
    pooled_std = math.sqrt(((n1 - 1) * group1.var(ddof=1) + (n2 - 1) * group2.var(ddof=1)) / (n1 + n2 - 2))
    return diff / pooled_std if pooled_std > 0 else 0.0


def _benjamini_hochberg(p_values: Sequence[float]):
    np = _require_dependency("numpy", _np)
    p = np.asarray(p_values)
    n = p.size
    order = np.argsort(p)
    ranked = np.empty_like(p)
    cumulative = 0.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        cumulative = min(cumulative + p[order[i]] * n / rank, 1.0)
        ranked[order[i]] = cumulative
    return ranked


def compute_paired_statistics(
    data: 'pd.DataFrame',
    donor_key: str,
    subtype_key: str,
    lam_label: str,
    trm_label: str,
    value_columns: Sequence[str],
) -> Dict[str, Tuple[float, float]]:
    pd = _require_dependency("pandas", _pd)
    stats = _require_dependency("scipy.stats", _stats)
    summary = {}
    for column in value_columns:
        pivot = data.pivot_table(index=donor_key, columns=subtype_key, values=column)
        lam = pivot[lam_label].dropna()
        trm = pivot[trm_label].dropna()
        common = lam.index.intersection(trm.index)
        if common.empty:
            raise ValueError(f"No paired data available for {column}")
        stat, p_value = stats.ttest_rel(lam.loc[common], trm.loc[common])
        summary[column] = (stat, p_value)
    return summary


def plot_scenith_pairs(
    data: 'pd.DataFrame',
    cfg: ScenithConfig,
    stats_summary: Dict[str, Tuple[float, float]],
    output_dir: Path,
) -> None:
    pd = _require_dependency("pandas", _pd)
    plt = _require_dependency("matplotlib", _plt)
    np = _require_dependency("numpy", _np)
    metrics = [cfg.basal_column, cfg.glycolysis_column]
    fig, axes = plt.subplots(1, len(metrics), figsize=(5 * len(metrics), 4))
    if len(metrics) == 1:
        axes = [axes]
    for ax, metric in zip(axes, metrics):
        pivot = data.pivot_table(index=cfg.donor_column, columns=cfg.subtype_column, values=metric)
        lam = pivot[cfg.lam_label]
        trm = pivot[cfg.trm_label]
        for donor in pivot.index:
            ax.plot([0, 1], [lam.loc[donor], trm.loc[donor]], marker="o", color="#636363", alpha=0.7)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([cfg.lam_label, cfg.trm_label])
        ax.set_ylabel(metric.replace("_", " ").title())
        stat, p_value = stats_summary.get(metric, (np.nan, np.nan))
        ax.set_title(f"{metric}\npaired t p={p_value:.3g}")
    fig.tight_layout()
    output_path = output_dir / cfg.output_prefix
    fig.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"Saved SCENITH paired plot to {output_path.with_suffix('.pdf')}")


def _write_placeholder(path: Path, message: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(message)


def _write_pdf_placeholder(path: Path, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"PDF placeholder: {title}\n")


def _write_tsv(path: Path, header: Sequence[str], rows: Sequence[Sequence[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")


def _resolve_output_path(base: Path, candidate: Optional[Path]) -> Path:
    if candidate is None:
        raise ValueError("Expected a path, received None")
    candidate = Path(candidate)
    return candidate if candidate.is_absolute() else base / candidate


