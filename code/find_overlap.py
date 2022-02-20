import pandas as pd
import numpy as np
import time
import sys

inputfile=sys.argv[1]#"change_s_e_node.tsv"
outputfile=sys.argv[2]#"longest_per_og.tsv"

print("Start Step2.5: Find the longest segment for each OG: ", time.asctime(time.localtime(time.time())))

### find the longest sequence for each og after format diamond output and change the start and the end
def longest_seq_per_og():
    df1 = pd.read_csv(inputfile, sep="\t",
                      names=["d_idx", "ContigLen", "OG", "PI", "Start", "End", "Similarity", "Seq", "COGFC", "Node"])

    ### find the longest protein seq for each OG
    SeqName_se_dic = dict()

    for idx in range(0, df1.index.max() + 1):
        seqname = str(df1.loc[idx, "Seq"])
        og = str(df1.loc[idx, "OG"])
        pi = float(df1.loc[idx, "PI"])
        s2 = int(df1.loc[idx, "Start"])
        e2 = int(df1.loc[idx, "End"])
        if (seqname, og, pi) not in SeqName_se_dic.keys():  # if seq and og do not appear before; create one
            SeqName_se_dic[(seqname, og, pi)] = [(s2, e2)]
        else:  # seq and og already exit
            # get the s1, e1 and see the overlap between the new line and existing protein length
            # the file is pre-sorted so only need to consider the last pair in the value
            s1 = SeqName_se_dic[(seqname, og, pi)][-1][0]
            e1 = SeqName_se_dic[(seqname, og, pi)][-1][1]
            if s2 <= e1:  # two protein overlap
                s = s1
                e = max(e1, e2)
                SeqName_se_dic[(seqname, og, pi)].pop()
                SeqName_se_dic[(seqname, og, pi)] = SeqName_se_dic.get((seqname, og, pi)) \
                                                    + [(s, e)]
            else:  # two protein do not overlap
                SeqName_se_dic[(seqname, og, pi)] = SeqName_se_dic.get((seqname, og, pi)) \
                                                    + [(s2, e2)]
    
    with open(outputfile, "w") as f:
        for k, v in SeqName_se_dic.items():
            for pair in v:
                f.write("{:>40}\t{:>10}\t{:8.5f}\t{:d}\t{:d}\n".format(str(k[0]), str(k[1]), float(k[2]), \
                                                                       int(pair[0]), int(pair[1])))
    print("End Step2.5: Find the longest segment for each OG: ", time.asctime(time.localtime(time.time())))

longest_seq_per_og()