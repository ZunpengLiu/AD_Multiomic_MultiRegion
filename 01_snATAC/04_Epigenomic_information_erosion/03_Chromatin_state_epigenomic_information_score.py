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

def Zcore_info_score(ratio_matrix):    
    z_score_matrix = z_score_normalize(ratio_matrix)
    score = np.sum(z_score_matrix.values * sign_matrix, axis=1)
    return score

# prepare the data
workdir="./05_ErosionScore"
os.chdir(workdir)

ChromatinStates=['EnhA1', 'EnhA2', 'EnhBiv', 'EnhG1', 'EnhG2', 'EnhWk', 'Het', 'Quies', 'ReprPC', 'ReprPCWk', 'TssA', 'TssBiv', 'TssFlnk', 'TssFlnkD','TssFlnkU', 'Tx', 'TxWk', 'ZNF.Rpts']

sign_matrix=[1, 1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1]

# calculate the erosion score for the chromatin states defined in different brain regions
for x in ["BSS00369","BSS00371","BSS01124","BSS01125","BSS01126","BSS01271,"BSS01272","BSS00077","BSS00078"]:
    raw, ratio_matrix = load_data(x+".Count_Fractions.txt",ChromatinStates)
    raw["Zcore_Score."+x]=Zcore_info_score(ratio_matrix)

raw["Score"]=raw[["Zcore_Score.BSS00369","Zcore_Score.BSS00371","Zcore_Score.BSS01124","Zcore_Score.BSS01125","Zcore_Score.BSS01126","Zcore_Score.BSS01271","Zcore_Score.BSS01272","Zcore_Score.BSS00077","Zcore_Score.BSS00078"]].mean(axis=1)

raw.to_csv("Score.txt",sep="\t",index=False)
