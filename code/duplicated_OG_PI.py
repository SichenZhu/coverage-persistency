# sorted_COGFC in diamond & nonduplicated COGFC with PI
# generate file that contains all the COGFC_PI that is corresponding to DIAMOND raw output

import pandas as pd
import numpy as np
import time
import sys

inputfile1=sys.argv[1]#"sorted_COGFC_diamond.tsv"
inputfile2=sys.argv[2]#"nondup_COGFunC_OG_PI_dnum.tsv"
outputfile=sys.argv[3]#"idx_COGFunC_OG_PI.tsv"

print("Start Step2.3 duplicated OG PI: ", time.asctime(time.localtime(time.time())))
idx_COGFC_df = pd.read_csv(inputfile1, sep="\t",\
                         names=["d_idx", "dCOGFC"])
COGFC_PI_df = pd.read_csv(inputfile2, sep="\t",\
                         names=["uCOGFC", "OG", "PI", "repeated_times"])

uOG_lst = list(COGFC_PI_df["OG"])
uPI_lst = list(COGFC_PI_df["PI"])
repeated_times = list(COGFC_PI_df["repeated_times"])

dOG_lst = []
dPI_lst = []
for u_idx in range(len(repeated_times)):
    times = int(repeated_times[u_idx])
    dOG_lst.extend([uOG_lst[u_idx]]*times)
    dPI_lst.extend([uPI_lst[u_idx]]*times)

#with open("record.tsv", "w") as f:
#    for idx in range(len(dOG_lst)):
#        f.write("{}\t{:8.5f}\n".format(dOG_lst[idx], dPI_lst[idx]))

idx_COGFC_df["OG"] = dOG_lst
idx_COGFC_df["PI"] = dPI_lst

idx_COGFC_df.to_csv(outputfile, sep="\t", index=False, header=False)
print("End Step2.3 duplicated OG PI: ", time.asctime(time.localtime(time.time())))
