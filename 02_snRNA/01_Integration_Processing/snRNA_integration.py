#!/usr/bin/env python3

def integrate(adata, 
              output=None, 
              prefix="", 
              batch=None,
              hvg=0, 
              use_harmony=False,
              max_iter_harmony=50,
              use_bbknn=False, 
              min_dist=0.5, 
              leiden="leiden_res1", 
              resolution=1.0, 
              plot=None, 
              compression=6,
              **kwargs):
    """
    Integration and clustering pipeline for snRNA-seq data using Scanpy, with optional Harmony or BBKNN.
    """

    print("Importing libraries")
    import scanpy as sc
    import pandas as pd
    import numpy as np
    import os
    import sys

    # Optional: set thread limits
    os.environ["OMP_NUM_THREADS"] = "30"
    os.environ["OPENBLAS_NUM_THREADS"] = "30"

    # Set figure resolution
    sc.settings.set_figure_params(dpi=350)

    # 1. Load data
    print("Loading data")
    adata = sc.read_h5ad(adata)
    print(f"Working with {adata.shape[0]} cells")

    # Store raw counts
    if "raw" in adata.layers:
        print('Copying .layers["raw"] to .X')
        adata.X = adata.layers["raw"].copy()
    else:
        print('Copying .X to .layers["raw"]')
        adata.layers["raw"] = adata.X.copy()

    # 2. Normalize and log-transform
    print("Normalizing and log-transforming data")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # 3. HVG selection
    if batch not in adata.obs.columns:
        batch = None
    elif len(pd.unique(adata.obs[batch])) == 0:
        batch = None
        print("No valid batch info found.")

    print("Identifying highly variable genes")
    if hvg > 0:
        sc.pp.highly_variable_genes(adata, n_top_genes=hvg, batch_key=batch, subset=True)
    else:
        sc.pp.highly_variable_genes(adata, min_mean=0.02, max_mean=4, min_disp=0.5)

    # 4. PCA
    print("Running PCA")
    sc.pp.pca(adata)
    rep = "X_pca"

    # 5. Batch correction
    if batch and use_harmony:
        print("Applying Harmony integration")
        sc.external.pp.harmony_integrate(adata, batch_key=batch, max_iter_harmony=max_iter_harmony)
        rep = "X_pca_harmony"

    if batch and use_bbknn:
        print("Applying BBKNN batch correction")
        sc.external.pp.bbknn(adata, batch_key=batch, use_rep=rep)
    else:
        print("Computing neighbors")
        sc.pp.neighbors(adata, use_rep=rep)

    # 6. UMAP
    print("Computing UMAP")
    sc.tl.umap(adata, min_dist=min_dist)

    # 7. Clustering
    print("Running Leiden clustering")
    sc.tl.leiden(adata, resolution=resolution, key_added=leiden)

    # 8. QC plots
    default_plots = [
        leiden, 'Region', 'Sample', 'batch', 'n_genes_by_counts', 
        'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 
        'total_counts_mt', 'log1p_total_counts_mt', 'pct_counts_mt', 
        'total_counts_ribo', 'log1p_total_counts_ribo'
    ]

    if plot is not None:
        plot = np.union1d(default_plots, plot)
    else:
        plot = default_plots

    print("Generating UMAP plots")
    for col in plot:
        if col in adata.obs.columns:
            print(f"Plotting {col}")
            sc.pl.umap(adata, color=col, save=f".{prefix}.{col}.png")
            sc.pl.umap(adata, color=col, legend_loc="on data", legend_fontsize=5, save=f".{prefix}.{col}.ondata.png")

    print("Generating violin plots")
    violin_metrics = ["pct_counts_mt", "pct_counts_ribo", "doublet_score", "log1p_total_counts"]
    for vp in np.intersect1d(violin_metrics, adata.obs.columns):
        sc.pl.violin(adata, vp, groupby=leiden, save=f".{prefix}.{leiden}_{vp}.png")

    # 9. Save output
    if output is not None:
        print("Saving results")
        adata.X = adata.layers["raw"].copy()
        del adata.layers["raw"]
        adata.write_h5ad(f"{output}{prefix}.h5ad", compression="gzip", compression_opts=compression)

    return adata


if __name__ == "__main__":
    from optparse import OptionParser

    usage = "\nWorkflow for snRNA-seq analysis using Scanpy (+ Harmony or BBKNN)\n"
    parser = OptionParser(usage, version="1.0")

    parser.add_option("-i", "--indir", dest="h5ad", type="string", help="Input .h5ad file")
    parser.add_option("-o", "--output", dest="output", type="string", help="Output directory prefix")
    parser.add_option("--prefix", default="", type=str, help="Prefix for saved files")
    parser.add_option("-b", "--batch", default=None, type=str, help="Batch column name")
    parser.add_option("--hvg", default=0, type=int, help="Number of HVGs to use (0 = auto)")
    parser.add_option("--use_harmony", action="store_true", default=False, help="Use Harmony for integration")
    parser.add_option("--max_iter_harmony", default=50, type=int, help="Max iterations for Harmony")
    parser.add_option("--use_bbknn", action="store_true", default=False, help="Use BBKNN for integration")
    parser.add_option("--min_dist", default=0.5, type=float, help="Minimum UMAP distance")
    parser.add_option("-l", "--leiden", default="leiden_clust", type=str, help="Leiden cluster column name")
    parser.add_option("-r", "--resolution", default=1.0, type=float, help="Leiden resolution")
    parser.add_option("--plot", default=None, type=str, help="Comma-separated list of obs keys to plot")
    parser.add_option("--compression", default=6, type=int, help="Compression level for output")

    (options, args) = parser.parse_args()

    integrate(
        options.h5ad,
        options.output,
        options.prefix,
        options.batch,
        options.hvg,
        options.use_harmony,
        options.max_iter_harmony,
        options.use_bbknn,
        options.min_dist,
        options.leiden,
        options.resolution,
        options.plot,
        options.compression
    )
