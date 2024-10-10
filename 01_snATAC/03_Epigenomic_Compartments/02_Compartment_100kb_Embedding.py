import pandas as pd
import scanpy as sc
import numpy as np
import os
import sys
import muon as mu
from muon import atac as ac
from matplotlib.pyplot import rc_context
from pandas.api.types import CategoricalDtype
import matplotlib.pyplot as plt
import seaborn as sns

# Set the working directory to the location where the data is stored
os.chdir("/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_ADMR_multiome_snRNA_snATAC/02_multiome_snATAC/06_HMM/100kb_emdeddings")

# Set visualization parameters for Scanpy
sc.figure_size = (8, 8)
sc.set_figure_params(scanpy=True, dpi=300, dpi_save=300, fontsize=6, color_map="viridis")

# Increase pandas display options to show more rows
pd.options.display.max_rows = 999

########## ----------------
# 1. Helper Functions
########## ----------------

# Function to plot UMAPs for a list of features
def pl_umap(adata, lsts, prefix):
    for lst in lsts:
        print(lst)
        sc.pl.umap(adata, color=lst, 
                   frameon=False, title="", 
                   legend_loc=None, 
                   vmin="p30", vmax="p95", 
                   show=False)
        plt.savefig(f"UMAP.{prefix}.{lst}.png", 
                    dpi=300, bbox_inches='tight', transparent=True)

# Function to generate violin plots for a list of features
def pl_violin(adata, lsts, prefix):
    for lst in lsts:
        sc.pl.violin(adata, groupby=lst, save=f".{prefix}.{lst}.png")


########## ----------------
# 2. Load the Data
########## ----------------

# Load the high-quality ATAC data in h5ad format
atac1 = sc.read("ATAC.Celltype_100kb_matrix.high_quality_100kb.h5ad")


########## ----------------
# 3. Transpose the Matrix
########## ----------------

# Transpose the matrix to switch from cell x 100kb to 100kb x cell format
atac1_T = atac1.copy().T

# Save the transposed matrix
atac1_T.write("ATAC.Celltype_100kb_matrix.high_quality_100kb.T.h5ad", compression="gzip")


########## ----------------
# 4. Run External Processing
########## ----------------

# External processing to be run in terminal (this step is not executed in the script)
'''
python snATAC_Integration_processing.py \
    -i ATAC.Celltype_100kb_matrix.high_quality_100kb.T.h5ad \
    -o ./ --prefix ATAC.Compartment_100kb.T. \
    --resolution 2
'''


########## ----------------
# 5. Transpose Matrix by Cell Type
########## ----------------

# Transpose and save the matrix separately for each major cell type
for x in atac1.obs["Cell_Class"].unique().tolist():
    print(x)
    tmp = atac1[atac1.obs["Cell_Class"].isin([x])]
    tmp.write(f"atac1.100kb_HighQualityBins.{x}.h5ad")

for x in atac1.obs["Cell_Class"].unique().tolist():
    print(x)
    tmp = atac1[atac1.obs["Cell_Class"].isin([x])]
    tmp2 = tmp.copy().T
    tmp2.write(f"atac1.100kb_HighQualityBins.{x}.T.h5ad")


########## ----------------
# 6. Add Epigenomic Signatures and Visualize
########## ----------------

# Load annotations for genomic compartments
annotations = pd.read_csv("../Compartment_Groups.with_Annotations.27993bins.tsv", sep="\t", header=0)
annotations.index = annotations["chr"] + ":" + (annotations["start"] - 1).astype(str) + "-" + annotations["end"].astype(str)

# Load the transposed matrix with clustering results
atac2 = sc.read("ATAC.Compartment_100kb.T.h5ad")

# Map the annotations to the ATAC dataset and visualize features
features_to_plot = ['LINE.L1', 'LINE.L2', 'SINE.Alu', 'SINE.MIR', 'LTR.ERV1', 'LTR.ERVK',
                    'LTR.ERVL', 'LTR.ERVL.MaLR', 'Satellite.acro', 'Satellite.centr',
                    'Satellite.telo', 'rRNA', 'tRNA', 'EnhA1', 'EnhA2', 'EnhBiv', 'EnhG1',
                    'EnhG2', 'EnhWk', 'Het', 'Quies', 'ReprPC', 'ReprPCWk', 'TssA',
                    'TssBiv', 'TssFlnk', 'TssFlnkD', 'TssFlnkU', 'Tx', 'TxWk', 'ZNF.Rpts',
                    'GW20_Cortex_LaminB', 'GW20_Cortex_SON', 'Compartment_Groups']

for x in features_to_plot:
    atac2.obs[x] = atac2.obs.index.map(annotations[x])

# Visualize the epigenomic signatures using UMAP
pl_umap(atac2, features_to_plot, "ATAC.Compartment_100kb.T")

# Visualize the 25 epigenomic compartment groups with custom colors
atac2.uns["Compartment_Groups_colors"] = ["#A00746", "#AE265F", "#BC4679", "#CB6593", "#D985AD", "#E7A4C7", "#F6C4E1", "#F5CCE4", 
                                          "#F5D4E7", "#F5DCEB", "#F5E4EE", "#F5ECF1", "#F5F5F5", "#D8E2F2", "#BBD0EF", "#9EBEED", 
                                          "#81ACEA", "#649AE7", "#4888E5", "#3F7CD8", "#3670CB", "#2D64BF", "#2458B2", "#1B4CA5", 
                                          "#134099"]

# Plot UMAP for the Compartment_Groups
pl_umap(atac2, ['Compartment_Groups'], "ATAC.Compartment_100kb.T")


########## ----------------
# 7. In Silico Correlation Analysis for Chromosomes
########## ----------------

# Define a function to compute correlation between chromosome features
def correlation_features(adata, feature):
    # Perform PCA and convert to dataframe
    df = pd.DataFrame(adata.obsm["X_lsi"])
    
    # Add annotations for the feature (e.g., chromosome) and group by it
    df["feature"] = adata.obs[feature].values
    df = df.groupby("feature").mean()
    
    # Compute correlation matrix
    corr = df.T.corr()
    return df, corr

# Perform correlation analysis on chromosome features
df = pd.DataFrame(atac2.obsm["X_lsi"])
feature = "Chromosome"
df["feature"] = atac2.obs[feature].values
df, corr = correlation_features(atac2, "Chromosome")

# Save the correlation matrix
corr.to_csv("atac2.Chromosome.correlation_matrix.tsv", sep="\t")

# Set custom categorical ordering for chromosomes (1-22 and X)
cat_order = CategoricalDtype(["chr" + str(x) for x in range(1, 23)] + ["chrX"], ordered=True)
atac2.obs["Chromosome"] = atac2.obs["Chromosome"].astype(cat_order)

# Define custom colors for the chromosomes and plot UMAP
atac2.uns["Chromosome_colors"] = ['#023fa5', '#4a6fe3', '#8595e1', '#7d87b9', '#b5bbe3', '#bec1d4',
                                  '#0fcfc0', '#9cded6', '#d5eae7', '#c6dec7', '#8dd593', '#11c638',
                                  '#ead3c6', '#f0b98d', '#ef9708', '#f3e1eb', '#f6c4e1', '#e6afb9',
                                  '#bb7784', '#e07b91', '#d33f6a', '#A00746', '#C6C6C6']

# Plot UMAP for Chromosome annotations
pl_umap(atac2, ['Chromosome'], "ATAC.Compartment_100kb.T")