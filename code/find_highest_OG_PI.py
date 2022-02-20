### find highest PI for each alignment
import pandas as pd
import numpy as np
import time
import sys

inputfile1=sys.argv[1]#"sorted_nondup_COGFunC_OG.tsv"
inputfile2=sys.argv[2]#"COG_sort_PI.tsv"
inputfile3=sys.argv[3]#"nCOG_PI.tsv"
outputfile=sys.argv[4]#"nondup_COGFunC_OG_PI.tsv"

print("Start Step2.3 find highest_OG_PI: ", time.asctime(time.localtime(time.time())))
COGFC_OG_df = pd.read_csv(inputfile1, sep="\t",\
                         names=["COGFC", "OG"])

COG_PI_df = pd.read_csv(inputfile2, sep="\t",\
                         names=["COG", "PI"])
nonOG_PI_df = pd.read_csv(inputfile3, sep="\t",\
                         names=["nCOG", "PI"])

PI_lst = []
newOG_lst = []

def find_highest_PI_for_COGFC(row):
    OG_lst = row["OG"].split(",")
    all_COG = True
    for OG in OG_lst: 
        if "COG" not in OG:
            all_COG = False
            break
    if all_COG:
        for idx in COG_PI_df.index:
            if COG_PI_df.loc[idx, "COG"] in OG_lst:
                PI_lst.append("{:7.5f}".format(COG_PI_df.loc[idx, "PI"]))
                newOG_lst.append(str(COG_PI_df.loc[idx, "COG"]))
                break
    else:
        lst_PI = []
        for OG in OG_lst:
            if "COG" in OG:
                lst_PI.append(float(COG_PI_df.loc[COG_PI_df["COG"]==OG, \
                                                  "PI"].to_numpy(\
                                                                 dtype=np.float, \
                                                                 copy=True)))
            else:
                lst_PI.append(float(nonOG_PI_df.loc[nonOG_PI_df["nCOG"]==OG,\
                                                   "PI"].to_numpy(\
                                                                 dtype=np.float,\
                                                                 copy=True)))
        idx = lst_PI.index(max(lst_PI))
        PI_lst.append("{:7.5f}".format(lst_PI[idx]))
        newOG_lst.append(str(OG_lst[idx]))

COGFC_OG_df.apply(find_highest_PI_for_COGFC, axis=1)

COGFC_OG_df["PI"] = PI_lst
COGFC_OG_df["OG"] = newOG_lst
COGFC_OG_df.to_csv(outputfile, sep="\t", index=False, header=False)

print("End Step2.3 find highest_OG_PI: ", time.asctime(time.localtime(time.time())))
