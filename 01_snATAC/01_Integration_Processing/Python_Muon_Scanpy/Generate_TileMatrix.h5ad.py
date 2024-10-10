#!/usr/bin/env python3

def generate_h5ad(indir, sample, peaks, outdir, **kwargs):
    import pandas as pd
    import numpy as np
    import anndata
    import muon as mu
    from muon import atac as ac
    import scanpy as sc
    import os
    import sys

    print("Reading fragments")
    os.chdir(indir)
    fragments=sample+".atac_fragments.tsv.gz"
    # read and ignore the header startswith #
    fragement_df = pd.read_csv(indir+fragments,sep="\t", header=None, index_col=None, engine="c",comment='#')

    fragement_df.columns = ["seqnames", "start", "end", "barcode", "count"]

    # groupby the fragments and sum the count
    fragement_df2=fragement_df.groupby("barcode")['count'].sum().reset_index()
    fragement_df2.columns = ["barcode","total"]
    print("Finished reading Fragments") 
    # cell metadata
    cm=fragement_df2
    cm["Sample"]=sample
    cm["Sample_barcode"]=cm["Sample"]+"#"+cm["barcode"]
    cm.index=cm["barcode"]
    cm.index.name = None
    cm.head()
    # annData
    atac = anndata.AnnData(obs=cm)

    # Read peak files
    print("Reading peak files")
    pk = pd.read_csv(peaks, sep="\t", header=0).iloc[:, [0,1,2]]
    pk.columns = ["seqnames", "start", "end"]
    pk.index = ["%s:%d-%d" % (chrom, begin, end) for chrom, begin, end in zip(pk["seqnames"], pk["start"], pk["end"])]
    pk["interval"] = pk.index.values.astype(str)

    # count fragments
    ac.tl.locate_fragments(atac,indir+fragments)

    print("Counting fragments")
    atac = ac.tl.count_fragments_features(atac, 
                                       pk.rename({"seqnames": "Chromosome", "start": "Start", "end": "End"}, axis=1), 
                                       extend_upstream=0, 
                                       extend_downstream=0)
    
    print("Converting to int16")
    atac.X = atac.X.astype(np.int16)
    # rename the obs index
    atac.obs.index=atac.obs["Sample_barcode"]
    cm.index.name = None
    # QC and Doublet calculation
    #print("QC of cells")
    #sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
    print(f"Before: {atac.n_obs} cells")
    mu.pp.filter_obs(atac, 'total', lambda x: (x >=1000))
    print(f"(After total_counts: {atac.n_obs} cells)")
        
    print("Writing h5ad file")
    atac.write_h5ad(outdir+sample+".h5ad", compression="gzip", compression_opts=9)

    print("Done!")


if __name__ == "__main__":
    from optparse import OptionParser

    MY_USAGE='\nWokrflow for snATAC-seq data analysis - Generating h5ad file (based on muon, scanpy, etc.)\n'
    parser=OptionParser(MY_USAGE,version='Version 1.0')
    parser.add_option('-i','--indir',   dest='indir',   type='string', help='Input directory')
    parser.add_option('-s','--sample',   dest='sample',   type='string', help='Sample name')
    parser.add_option('-p','--peak',   dest='peaks',   
                      type='string', help='Peak file',
                      default="~/hg38_5000bp.bed")
    parser.add_option('-o','--outdir',  dest='outdir',  type='string', help='Output directory')
    generate_h5ad(**vars(parser.parse_args()[0]))
