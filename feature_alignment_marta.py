
import os
from ms_preprocessing_oms import fill_feature_map, feature_alignment_oms
import pandas as pd

paths = [
r"D:\file1.csv",
r"D:\file2.csv"
]

sample_names = [os.path.basename(p).split('.')[0] for p in paths]
dfs = [pd.read_csv(p) for p in paths]

fms = []
for n in range(len(dfs)):
    fms.append(fill_feature_map(dfs[n]['mz'].values, 
                                dfs[n]['rt'].values, 
                                dfs[n]['intensity'].values))

df_align, cons_map = feature_alignment_oms(fms, sample_names)

print(df_align)