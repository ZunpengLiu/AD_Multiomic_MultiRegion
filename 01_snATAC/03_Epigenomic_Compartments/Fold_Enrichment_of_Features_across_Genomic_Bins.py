import pandas as pd
import numpy as np
import os
import sys

# List all files ending with "_18_CALLS_segments.bed"
files = os.listdir("./")
files = [i for i in files if i.endswith("_18_CALLS_segments.bed")]

# Process each segment file
for i in files:
    # Intersect the hg38 100kb high-quality regions with the current segment file
    os.system(f"bedtools intersect -a hg38_100kb.HighQuality.bed -b {i} -wao > hg38_100kb.{i.split('_18_CALLS_segments.bed')[0]}.count.txt")
    
    name = i.split("_18_CALLS_segments.bed")[0]
    count = pd.read_csv(f"hg38_100kb.{name}.count.txt", sep="\t", header=None)
    
    # Group by columns ["0", "1", "2", "6"] and sum the counts in column [12]
    count = count.groupby([0, 1, 2, 6])[12].sum().reset_index()
    count.columns = ["chr", "start", "end", "state", "count"]

    # Calculate state percentages
    state = pd.read_csv(i, sep="\t", header=None)
    state["length"] = state[2] - state[1] + 1  # Add a length column for each region
    pct = state.groupby(3)["length"].sum().reset_index()
    pct["pct"] = pct["length"] / pct["length"].sum()  # Calculate percentage
    pct.set_index(3, inplace=True)

    # Map the state percentages and calculate overlap percentage
    count["state_pct"] = count["state"].map(pct["pct"])
    count["overlap_pct"] = count["count"] / 100000  # Normalize by region size (100kb)
    
    # Calculate fold enrichment and log2 fold enrichment
    count["fold_enrichment"] = count["overlap_pct"] / count["state_pct"]
    count["log2_fold_enrichment"] = np.log2(count["fold_enrichment"])

    # Reshape to wide format (one row per region, one column per state)
    count = count.pivot_table(index=["chr", "start", "end"], columns="state", values="log2_fold_enrichment").reset_index()
    count.index = count["chr"] + ":" + count["start"].astype(str) + "-" + count["end"].astype(str)
    count.drop(["chr", "start", "end"], axis=1, inplace=True)
    count.fillna(0, inplace=True)  # Replace any missing values with 0

    # Save log2 fold enrichment data
    count.to_csv(f"hg38_100kb.{name}.log2_fold_enrichment.tsv", sep="\t")

# Combine enrichment data for different brain regions
df_mean = pd.read_csv("hg38_100kb.BSS01125.log2_fold_enrichment.tsv", sep="\t", index_col=0)
df_mean.columns = ["BSS01125_" + j for j in df_mean.columns]

# Combine data across multiple samples
for i in ["BSS01126", "BSS00369", "BSS00371", "BSS01271", "BSS01272", "BSS00077", "BSS00078"]:
    tmp = pd.read_csv(f"hg38_100kb.{i}.log2_fold_enrichment.tsv", sep="\t", index_col=0)
    tmp.columns = [i + "_" + j for j in tmp.columns]
    df_mean = df_mean.merge(tmp, left_index=True, right_index=True, how="inner")

# Calculate the mean fold enrichment for specific states across samples
states = ["EnhA1", "EnhA2", "EnhBiv", "EnhG1", "EnhG2", "EnhWk", "Het", "Quies", "ReprPC", "ReprPCWk", 
          "TssA", "TssBiv", "TssFlnk", "TssFlnkD", "TssFlnkU", "Tx", "TxWk", "ZNF/Rpts"]

for x in states:
    df_mean[x] = df_mean[[i + "_" + x for i in ["BSS01125", "BSS01126", "BSS00369", "BSS00371", "BSS01271", "BSS01272", "BSS00077", "BSS00078"]]].mean(axis=1)

# Save the combined log2 fold enrichment data
df_mean2 = df_mean[states]
df_mean2.to_csv("hg38_100kb.chromHMM_mean.log2_fold_enrichment.tsv", sep="\t")

############################## Lamin B1 and SON Analysis

# List all relevant files for Lamin B1 and SON analysis
files = os.listdir("/net/bmc-lab5/data/kellis/group/Zunpeng_Sharing/LAD_SON")
files = [i for i in files if i.endswith("_peaks_intersect_srt.bed")]

# Process each file
for i in files:
    os.system(f"bedtools intersect -a hg38_100kb.HighQuality.bed -b /net/bmc-lab5/data/kellis/group/Zunpeng_Sharing/LAD_SON/{i} -wao > hg38_100kb.{i.split('_peaks_intersect_srt.bed')[0]}.count.txt")
    
    name = i.split("_peaks_intersect_srt.bed")[0]
    count = pd.read_csv(f"hg38_100kb.{name}.count.txt", sep="\t", header=None)
    count = count.groupby([0, 1, 2])[9].sum().reset_index()
    count.columns = ["chr", "start", "end", "count"]
    count["state"] = name

    state = pd.read_csv(f"/net/bmc-lab5/data/kellis/group/Zunpeng_Sharing/LAD_SON/{i}", sep="\t", header=None)
    state["length"] = state[2] - state[1] + 1
    state["state"] = name
    pct = state.groupby("state")["length"].sum().reset_index()
    pct["pct"] = pct["length"] / 3137300923  # Percentage of total genome length (3.13 billion bp)
    pct.set_index("state", inplace=True)

    count["state_pct"] = count["state"].map(pct["pct"])
    count["overlap_pct"] = count["count"] / 100000
    count["fold_enrichment"] = count["overlap_pct"] / count["state_pct"]
    count["log2_fold_enrichment"] = np.log2(count["fold_enrichment"] + 1e-10)  # Add a small value to avoid log(0)

    # Reshape to wide format and save
    count = count.pivot_table(index=["chr", "start", "end"], columns="state", values="log2_fold_enrichment").reset_index()
    count.index = count["chr"] + ":" + count["start"].astype(str) + "-" + count["end"].astype(str)
    count.drop(["chr", "start", "end"], axis=1, inplace=True)
    count.fillna(0, inplace=True)
    count.to_csv(f"hg38_100kb.{name}.log2_fold_enrichment.tsv", sep="\t")

# Combine Lamin B1, SON, and ChromHMM data

# Directory containing enrichment data
indir = '~/Cortex_LaminB1_SON_H3K9me2'
os.chdir(indir)

# Read in all fold enrichment files
files = [x for x in os.listdir(indir) if x.endswith(".log2_fold_enrichment.tsv")]

# Initialize an empty DataFrame
df = pd.DataFrame()

# Process and concatenate all fold enrichment files
for i in files:
    print(f"Processing {i}")
    name = i.replace(".log2_fold_enrichment.tsv", "").replace("hg38_100kb.", "")
    tmp = pd.read_csv(i, sep="\t", header=None)
    tmp.columns = ["chr.start.end", "value"]
    tmp["celltype"] = name
    df = pd.concat([df, tmp])

# Pivot to wide format
df2 = df.pivot_table(index=["chr.start.end"], columns="celltype", values="value").reset_index()

# Save the combined data
df2.to_csv("LAD_SON_H3K9me2.100kb.value.tsv.gz", sep="\t", compression="gzip")

# Final combination of repeats, ChromHMM, LaminB, and SON data

# Read in different datasets
repeats = pd.read_csv("hg38_100kb.Repeats.log2_fold_enrichment.tsv", sep="\t", index_col=0)
chromhmm = pd.read_csv("hg38_100kb.chromHMM_mean.log2_fold_enrichment.tsv", sep="\t", index_col=0)
lads = pd.read_csv("hg38_100kb.GW20_Cortex_LaminB.log2_fold_enrichment.tsv", sep="\t", index_col=0)
Son = pd.read_csv("hg38_100kb.GW20_Cortex_SON.log2_fold_enrichment.tsv", sep="\t", index_col=0)

# Combine all datasets
combine = repeats.merge(chromhmm, left_index=True, right_index=True, how="inner")
combine = combine.merge(lads, left_index=True, right_index=True, how="inner")
combine = combine.merge(Son, left_index=True, right_index=True, how="inner")

# Save the combined dataset
combine.to_csv("hg38_100kb.Repeats.ChromHMM_LaminB_SON.log2_fold_enrichment_all.tsv", sep="\t")
