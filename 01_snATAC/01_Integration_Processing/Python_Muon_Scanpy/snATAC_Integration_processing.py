#!/usr/bin/env python3
# Integration of h5ad files

def integrate(adata, 
              output=None, 
              prefix="",
              batch=None,
              min_n_cells_by_counts=20, 
              hvg=0,
              use_harmony=False,
              max_iter_harmony=50,
              use_bbknn=False,
              resolution=1.0, 
              min_dist=0.5, 
              leiden="leiden", 
              plot=None, 
              compression=6,
              **kwargs):
    
    
    import pandas as pd
    import numpy as np
    import muon as mu
    from muon import atac as ac
    import scanpy as sc
    import os
    import sys
    #import bbknn

    sc.settings.set_figure_params(dpi=350)

    atac=sc.read_h5ad(adata)

    if batch not in atac.obs.columns:
        batch = None
    
    print("Remove samples less than 5 cells detected for the given batch")
    count_of_cells = atac.obs.groupby(batch).size()
    atac = atac[atac.obs[batch].isin(count_of_cells[count_of_cells >= 5].index.tolist())]
    
    if "raw" not in atac.layers:
        print("Copying .X to .layers[\"raw\"]")
        atac.layers["raw"] = atac.X.copy()
    else:
        print("Copying .layers[\"raw\"] to .X")
        atac.X = atac.layers["raw"].copy()

    print("Re-calculating qc metrics")
    #sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=True, inplace=True, layer="raw")

    if min_n_cells_by_counts > 0 and "n_cells_by_counts" in atac.var.columns:
        print("Filter peaks which ATAC signal is not well detected")
        mu.pp.filter_var(atac, "n_cells_by_counts", lambda x: x >= min_n_cells_by_counts)

    #TF-IDF normalisation is implemented in the muon’s ATAC module
    print("Running TF-IDF normalisation")
    ac.pp.tfidf(atac, scale_factor=1e4)

    print("Running LSI to get latent components with TF-IDF counts")
    ac.tl.lsi(atac)
    atac.X = atac.layers["raw"].copy()
    atac.X = atac.X.astype(np.float32)
    del atac.layers["raw"]

    print("Removing first component (which is typically associated with number of peaks or counts per cell)")
    # We find the first component is typically associated with number of peaks or counts per cell so it is reasonable to remove it:
    atac.obsm['X_lsi'] = atac.obsm['X_lsi'][:,1:]
    atac.varm["LSI"] = atac.varm["LSI"][:,1:]
    atac.uns["lsi"]["stdev"] = atac.uns["lsi"]["stdev"][1:]
    

    use_rep="X_lsi"
    n_pcs=len(atac.uns["lsi"]["stdev"])
    #atac.write_h5ad(output+prefix+".tmp.h5ad", compression="gzip", compression_opts=compression)
    if batch is not None and use_harmony:
        print("Running Harmony")
        use_rep_adj = "%s_harmony" % use_rep
        sc.external.pp.harmony_integrate(atac, batch, 
                                         basis=use_rep, 
                                         adjusted_basis=use_rep_adj, 
                                         max_iter_harmony=max_iter_harmony)
        use_rep = use_rep_adj
    
    if batch is not None and use_bbknn:
        print("Running BBKNN")
        sc.external.pp.bbknn(atac, batch, use_rep=use_rep, n_pcs=n_pcs)
        #bbknn.bbknn(atac, batch, use_rep=use_rep, n_pcs=n_pcs)

    else:
        print("Running sc.pp.neighbors nearest neighbors with use_rep: "+use_rep)
        sc.pp.neighbors(atac, use_rep=use_rep, n_pcs=n_pcs)

    print("Computing UMAP")
    sc.tl.umap(atac, min_dist=min_dist)

    print("Clustering cells by leiden")
    sc.tl.leiden(atac, key_added=leiden, resolution=resolution)

    # 10. Run marker genes and dotplot
    print("Plotting")
    default_plots=[ "Chr",'TSSEnrichment', 'ReadsInTSS', 'ReadsInPromoter', 
       'ReadsInBlacklist', 'PromoterRatio', 'PassQC', 'NucleosomeRatio', 
       'nMultiFrags', 'nMonoFrags', 'nFrags', 'nDiFrags', 'DoubletScore', 
       'DoubletEnrichment', 'BlacklistRatio', 'n_genes_by_counts', 
       'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 
       'total_counts_Distal', 'log1p_total_counts_Distal', 'pct_counts_Distal', 
       'total_counts_Promoter', 'log1p_total_counts_Promoter', 'pct_counts_Promoter', 
       'total_counts_Exonic', 'log1p_total_counts_Exonic', 'pct_counts_Exonic', 
         'region', 'frag']

    
    if plot is not None:
        plot = np.union1d(default_plots, plot)
    else:
        plot = default_plots

    # set up the resolution dpi for the plot
    sc.settings.set_figure_params(dpi=350)

    print("UMAP plots of basic QC metrics")
    for col in plot:
        if col in atac.obs.columns:
            print("Plotting %s" % col)
            sc.pl.umap(atac, color=col, save="."+prefix+".%s.png" % col)
            sc.pl.umap(atac, color=col, legend_loc="on data", legend_fontsize=5, save="."+prefix+".%s.ondata.png" % col)

    print("Plotting %s" % leiden)
    sc.pl.umap(atac, color=leiden, save="_%s_beside.png" % leiden)
    sc.pl.umap(atac, color=leiden, save="_%s_ondata.png" % leiden, legend_loc="on data", legend_fontsize=4)

    print("Writing to h5ad")
    #atac.write_h5ad(output+prefix+".h5ad", compression="gzip", compression_opts=compression)
    atac.write_h5ad(output+prefix+".h5ad", compression="gzip")



if __name__ == "__main__":
    from optparse import OptionParser

    MY_USAGE='\nWokrflow for snATAC-seq data analysis - Generating h5ad file (based on snapATAC2, muon, scanpy)'
    parser=OptionParser(MY_USAGE,version='Version 1.0')
    parser.add_option('-i','--input', dest='h5ad',   type='string', help='Input h5ad file')
    parser.add_option("-o", "--output", dest="output", type="string", help="Output directory of the h5ad file")
    parser.add_option('-p','--prefix',   dest='prefix',   type='string', default="", 
                      help='Prefix of leiden clustering result, which is used in adata.obs[leiden] = ["%s%s" % (prefix, v) for v in adata.obs[leiden].values.astype(str)]')
    parser.add_option('-b',"--batch", default=None, type=str, help='Batch column name')
    parser.add_option('--min_n_cells_by_counts', dest='min_n_cells_by_counts', type='int', default=20,
                        help='Minimum number of cells per gene by counts, which is used in sc.pp.filter_genes(adata, min_cells=min_n_cell_by_counts)')
    parser.add_option("--hvg", default=0, type=int)
    parser.add_option("--use_harmony", dest= "use_harmony", default=False, action="store_true")
    parser.add_option("--max_iter_harmony", default=50, type=int)    
    parser.add_option("--use_bbknn", dest= "use_bbknn", default=False, action="store_true")
    parser.add_option("-r", "--resolution", default=1.0, type=float, help="Resolution for Leiden")
    parser.add_option("--min_dist", default=0.5, type=float, help="Minimum distance for UMAP")
    parser.add_option("-l","--leiden", default="leiden_clust", type=str,help="Leiden clustering column")
    parser.add_option("--plot", default=None, type=str)
    parser.add_option("--compression", default=6, type=int)

    (options,args)=parser.parse_args()

    integrate(options.h5ad, options.output, options.prefix, options.batch, options.min_n_cells_by_counts, options.hvg,
              options.use_harmony, options.max_iter_harmony, 
              options.use_bbknn, options.resolution, options.min_dist, options.leiden, 
              options.plot, options.compression)
