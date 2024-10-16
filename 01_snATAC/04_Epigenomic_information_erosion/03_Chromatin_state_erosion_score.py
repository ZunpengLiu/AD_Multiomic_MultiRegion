import numpy as np
import pandas as pd
import os
import sys

pd.set_option('display.max_columns', 1000000)
pd.set_option('display.max_rows', 1000000)

# Function to load the data
def load_data(filename,ChromatinStates):
    # prepare the data
    print("loading data from ",filename)
    res=pd.read_csv(filename,sep="\t")
    raw=res
    res=res.loc[:,res.columns.str.contains(filename.split(".")[0])]
    res=res.loc[:,res.columns.str.contains("Ratio")]
    res.columns=res.columns.str.replace("_Ratio","")
    res.columns=res.columns.str.replace(filename.split(".")[0]+".","")
    ratio_matrix=res
    ratio_matrix=ratio_matrix.fillna(0)
    ratio_matrix=ratio_matrix.loc[:,ChromatinStates]
    return raw, ratio_matrix

# Function to normalize the data
def z_score_normalize(matrix):
    # Calculating the mean of each column
    column_means = np.mean(matrix, axis=0)
    
    # Calculating the standard deviation of each column
    column_stddevs = np.std(matrix, axis=0)
    
    # calculating the Z-score normalization for each element
    z_score_matrix = (matrix - column_means) / column_stddevs
    
    return z_score_matrix

def Zcore_erosion_score(ratio_matrix):    
    z_score_matrix = z_score_normalize(ratio_matrix)
    erosion_score = np.sum(z_score_matrix.values * sign_matrix, axis=1)
    return erosion_score

# prepare the data
workdir="./05_ErosionScore"
os.chdir(workdir)

ChromatinStates=['EnhA1', 'EnhA2', 'EnhBiv', 'EnhG1', 'EnhG2', 'EnhWk', 'Het', 'Quies', 'ReprPC', 'ReprPCWk', 'TssA', 'TssBiv', 'TssFlnk', 'TssFlnkD','TssFlnkU', 'Tx', 'TxWk', 'ZNF.Rpts']

sign_matrix=[-1, -1, -1, -1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, 1]

# calculate the erosion score for the chromatin states defined in different brain regions
for x in ["BSS00369","BSS01126","BSS01272","BSS00077"]:
    raw, ratio_matrix = load_data(x+".Count_Fractions.txt",ChromatinStates)
    raw["Zcore_ErosionScore."+x]=Zcore_erosion_score(ratio_matrix)

# BSS00369 for the frontal cortex, BSS01126 for the hippocampus, BSS01272 for the middle frontal area, and BSS00077 for the angular gyrus
raw["ErosionScore"]=raw[["Zcore_ErosionScore.BSS00369",
                            "Zcore_ErosionScore.BSS01126",
                            "Zcore_ErosionScore.BSS01272",
                            "Zcore_ErosionScore.BSS00077"]].mean(axis=1)

raw.to_csv("ErosionScore.txt",sep="\t",index=False)