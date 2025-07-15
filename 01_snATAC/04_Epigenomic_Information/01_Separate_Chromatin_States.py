import pandas as pd
import numpy as np
import os
import sys

files=os.listdir("./")
files=[i for i in files if i.endswith("_18_CALLS_segments.bed")]

for file in files:
    print("processing "+file)
    name=file.split("_18_CALLS_segments.bed")[0]
    df=pd.read_csv(file,sep="\t",header=None)

    for state in df[3].unique().tolist():
        tmp=df[df[3]==state]
        state=state.replace("/",".")
        tmp.to_csv(name+"/"+state+".bed",sep="\t",index=False,header=False)

