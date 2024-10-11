import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
import itertools
import math
import scipy
from sklearn.preprocessing import MinMaxScaler
from scipy.stats import entropy
import seaborn as sns
import matplotlib.pyplot as plt

from pandas.api.types import CategoricalDtype


pd.options.display.max_rows = 80000

######## Functions ########
# the following functions are adapted from the papers with some modifications:
# Li, Yang Eric, Sebastian Preissl, Michael Miller, Nicholas D. Johnson, Zihan Wang, Henry Jiao, Chenxu Zhu, et al. 2023. “A Comparative Atlas of Single-Cell Chromatin Accessibility in the Human Brain.” Science (New York, N.Y.) 382 (6667): eadf7044.
# Zu, Songpeng, Yang Eric Li, Kangli Wang, Ethan J. Armand, Sainath Mamde, Maria Luisa Amaral, Yuelai Wang, et al. 2023. “Single-Cell Analysis of Chromatin Accessibility in the Adult Mouse Brain.” Nature 624 (7991): 378–89.
def saveH(prefix, X):
    print("=== write matrix H ===")
    fileN = [prefix, "H", "mx"]
    fileN = ".".join(fileN)
    np.savetxt(fileN, X, fmt="%g", delimiter="\t")


def saveW(prefix, X):
    print("=== write matrix W ===")
    fileN = [prefix, "W", "mx"]
    fileN = ".".join(fileN)
    np.savetxt(fileN, X, fmt="%g", delimiter="\t")


def cal_featureScore_kim(W):
    """extract feature from W"""
    print("=== extract feature from W ===")
    k = W.shape[1]
    m = W.shape[0]
    s_list = []
    for i in range(m):
        rowsum = np.sum(W[i,])
        p_iq_x_list = []
        for q in range(k):
            p_iq = W[i, q] / rowsum
            if p_iq != 0:
                tmp = p_iq * math.log(p_iq, 2)
            else:
                tmp = 0
            p_iq_x_list.append(tmp)
        s = 1 + 1 / math.log(k, 2) * np.sum(p_iq_x_list)
        s_list.append(s)
    return s_list


def predict_H(H):
    """extract feature from H"""
    print("=== extract feature from H ===")
    colmax = np.amax(H, axis=0)
    colsum = np.sum(H, axis=0)
    p = colmax / colsum
    idx = H.argmax(axis=0)
    out = [idx, p]
    return out

def predict_H2(H):
    """extract feature from H"""
    print("=== extract feature from H ===")
    # Normalize H to scale [0, 1] by column
    scaler = MinMaxScaler()
    H_scaled = pd.DataFrame(scaler.fit_transform(H.T).T, 
                            columns=data.columns)
    H_scaled.index=["Module-"+str(i) for i in range(1,151)]
    # wide to long format with index_name
    H_scaled2 = H_scaled.melt(var_name='celltype', 
                              value_name='value', 
                              ignore_index=False).reset_index()
    H_scaled2["value"].describe()
    # get 95% quantile
    q = H_scaled2["value"].quantile(0.95)
    # Apply threshold to associate clusters with modules
    H_scaled3 = H_scaled2[H_scaled2["value"]>q]
    H_scaled3["Subtype"]=H_scaled3["celltype"].str.split(".").str[0]
    H_scaled3["BrainRegion"]=H_scaled3["celltype"].str.split(".").str[1]
    H_scaled3["Subtype"]
    H_wide=H_scaled3.pivot(index='index', columns='celltype', values='value')

    # rank modules by the counts in all celltype (columns)
    H_wide["counts"] = H_wide.count(axis=1)
    H_wide2=H_wide.sort_values(by="counts",ascending=False)

    return H_scaled3


####### Main #######
# load data
outdir="./04_Modules/01_Subtype_peaks"
os.chdir(outdir)
data = pd.read_csv("ATAC.Subtype_BrainRegion.PeakScoreMatrix.tsv",sep="\t",
                   index_col=0)

# Run NMF
prefix="Subtype_BrainRegion_M150"
n_components = 150
model = NMF(n_components=n_components, init='nndsvda', random_state=0)
W = model.fit_transform(data) 
H = model.components_ 

# save raw H and W
saveH(prefix, H)
saveW(prefix, W)