import pandas as pd
import numpy as np
import time
import sys

inputfile1=sys.argv[1]#"reformat_diamond_output.tsv"
outputfile=sys.argv[2]#"change_s_e.tsv"

print("Start Step2.4 change start end point in DIAMOND output: ", time.asctime(time.localtime(time.time())))
idx_COGFC_df = pd.read_csv(inputfile1, sep="\t",\
                         names=["d_idx", "ContigLen", \
                               "OG", "PI", \
                               "Start", "End", "Similarity", \
                               "Seq", "COGFC"])
#uOG_lst = idx_COGFC_df["OG"].unique()

for idx in range(0, idx_COGFC_df.index.max()+1):
    if idx_COGFC_df.loc[idx, "Start"] > idx_COGFC_df.loc[idx, "End"]:
        idx_COGFC_df.loc[idx, "Start"], idx_COGFC_df.loc[idx, "End"] = \
        idx_COGFC_df.loc[idx, "End"], idx_COGFC_df.loc[idx, "Start"]
idx_COGFC_df.to_csv(outputfile,sep="\t", index=False, header=False)
print("End Step2.4 change start end point in DIAMOND output: ", time.asctime(time.localtime(time.time())))
