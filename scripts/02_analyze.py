#!/usr/bin/env python3
"""
Single-cell analysis of GSE231935 (M3-Seq) E. coli data.
Pipeline: load -> QC -> normalize -> reduce -> cluster -> DE -> plots.
"""

import warnings
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from anndata import AnnData

warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "output"

# M3-Seq E. coli dataset (GSM7306272, scifi8, MG1655)
EXPR_FILE = DATA_DIR / "GSM7306272_expression_filtered_scifi8_MG1655_25.csv.gz"
CELL_INDEX = DATA_DIR / "GSM7306272_cell_index_scifi8_MG1655_25.csv.gz"
GENE_INDEX = DATA_DIR / "GSM7306272_gene_index_scifi8_MG1655_25.csv.gz"


def load_m3seq_ecoli() -> AnnData:
    """Load M3-Seq E. coli expression matrix into AnnData."""
    print("Loading expression matrix...")
    expr = pd.read_csv(EXPR_FILE, index_col=0)

    # Keep only E. coli genes (EC_)
    ec_genes = [c for c in expr.columns if c.startswith("EC_")]
    expr = expr[ec_genes]

    print("Loading cell metadata...")
    cell_meta = pd.read_csv(CELL_INDEX)
    if expr.shape[0] != len(cell_meta):
        raise ValueError(
            f"Cell count mismatch: expression {expr.shape[0]} vs cell_index {len(cell_meta)}"
        )

    print("Loading gene metadata...")
    gene_meta = pd.read_csv(GENE_INDEX)
    ec_gene_meta = (
        gene_meta[gene_meta["Name"].str.strip().str.startswith("EC_", na=False)]
        .set_index("Name")
        .reindex(expr.columns)
    )

    adata = AnnData(
        expr.values,
        obs=cell_meta.reset_index(drop=True),
        var=ec_gene_meta.fillna({"chr": "MG1655"}),
    )
    adata.obs_names = [f"cell_{i}" for i in range(adata.n_obs)]
    adata.var_names = expr.columns.tolist()

    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def run_analysis(adata: AnnData) -> AnnData:
    """Run Scanpy pipeline: QC, normalize, reduce, cluster, DE."""
    # QC
    def _flatten(x):
        return np.asarray(x).flatten()

    adata.var["n_cells"] = _flatten((adata.X > 0).sum(axis=0))
    adata.obs["n_genes"] = _flatten((adata.X > 0).sum(axis=1))
    adata.obs["total_counts"] = _flatten(adata.X.sum(axis=1))

    # Filter very low-quality cells (bacterial data is sparse)
    sc.pp.filter_cells(adata, min_genes=50)
    sc.pp.filter_genes(adata, min_cells=10)

    print(f"After QC: {adata.n_obs} cells x {adata.n_vars} genes")

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG (bacterial data may have few variable genes)
    sc.pp.highly_variable_genes(adata, n_top_genes=min(500, adata.n_vars))

    # Scale and reduce
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=min(50, adata.n_obs - 1, adata.n_vars - 1))

    # UMAP (use PCA for very sparse data)
    try:
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(30, adata.obsm["X_pca"].shape[1]))
        sc.tl.umap(adata)
    except Exception:
        print("UMAP failed (sparse data); using PCA for visualization")
        adata.obsm["X_umap"] = adata.obsm["X_pca"][:, :2]

    # Cluster (use KMeans on PCA if leiden/louvain unavailable)
    try:
        sc.tl.leiden(adata, resolution=0.5)
    except (ImportError, Exception):
        try:
            sc.tl.louvain(adata, resolution=0.5)
            adata.obs["leiden"] = adata.obs["louvain"]
        except (ImportError, Exception):
            from sklearn.cluster import KMeans
            n_clusters = min(8, adata.n_obs // 100)
            adata.obs["leiden"] = (
                KMeans(n_clusters=n_clusters, random_state=0, n_init=10)
                .fit_predict(adata.obsm["X_pca"][:, :min(30, adata.obsm["X_pca"].shape[1])])
                .astype(str)
            )

    # Differential expression
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")

    return adata


def plot_results(adata: AnnData) -> None:
    """Generate UMAP, violin, and heatmap figures."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(OUTPUT_DIR)

    # UMAP by cluster
    fig, ax = plt.subplots(figsize=(6, 5))
    if "X_umap" in adata.obsm:
        sc.pl.umap(adata, color="leiden", ax=ax, show=False, title="E. coli clusters (Leiden)")
    else:
        sc.pl.pca(adata, color="leiden", ax=ax, show=False, title="E. coli clusters (PCA)")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "umap_clusters.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {OUTPUT_DIR / 'umap_clusters.png'}")

    # Violin: total counts by cluster
    fig, ax = plt.subplots(figsize=(8, 4))
    sc.pl.violin(adata, keys="total_counts", groupby="leiden", ax=ax, show=False)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "violin_total_counts.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {OUTPUT_DIR / 'violin_total_counts.png'}")

    # Top DE genes heatmap (if we have enough clusters)
    n_clusters = adata.obs["leiden"].nunique()
    if n_clusters >= 2:
        try:
            sc.pl.rank_genes_groups_heatmap(
                adata,
                n_genes=5,
                groupby="leiden",
                show=False,
                save="_top5.png",
                figsize=(8, 6),
            )
            print(f"Saved heatmap to {OUTPUT_DIR}")
        except Exception as e:
            print(f"Heatmap skipped: {e}")

    # Dotplot of top DE genes per cluster
    if n_clusters >= 2:
        try:
            sc.pl.rank_genes_groups_dotplot(
                adata,
                n_genes=3,
                groupby="leiden",
                show=False,
                save="_top3.png",
                figsize=(6, 4),
            )
            print(f"Saved dotplot to {OUTPUT_DIR}")
        except Exception as e:
            print(f"Dotplot skipped: {e}")

    # Summary stats
    stats = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "n_clusters": adata.obs["leiden"].nunique(),
        "cluster_sizes": adata.obs["leiden"].value_counts().to_dict(),
    }
    pd.DataFrame([stats]).to_csv(OUTPUT_DIR / "summary_stats.csv", index=False)
    print(f"Saved {OUTPUT_DIR / 'summary_stats.csv'}")


def main():
    if not EXPR_FILE.exists():
        raise FileNotFoundError(
            f"Data not found. Run 01_download_data.py first.\nExpected: {EXPR_FILE}"
        )

    adata = load_m3seq_ecoli()
    adata = run_analysis(adata)
    plot_results(adata)
    print("Analysis complete.")


if __name__ == "__main__":
    main()
