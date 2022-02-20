import pandas as pd
import numpy as np
import time
import sys

inputfile=sys.argv[1]#"sorted_longest_per_og_len_node.tsv"
outputfile1=sys.argv[2]#"overlap.tsv"
outputfile2=sys.argv[3]#"nonoverlap.tsv"

print("Start Step2.6: find overlap & nonoverlap part between each OG: ", time.asctime(time.localtime(time.time())))

### find overlap between two different OGs (duplicated)
def overlap_btw_diff_OGs_unmodified():
    df1 = pd.read_csv(inputfile, sep="\t", \
                      names=["Seq", "OG", "PI", "Start", "End", "Len", "Node"])
    
    def Search(L, a, s, e):  # s and e is the index of the line_index in seq_linenum_dic.values()
        # the start point of each alignment is sorted
        def bSearch(L, a, low, high):
            if L[low] > a:
                return None
            elif L[high] <= a:
                return high
            else:
                mid = (low + high) // 2
                if high - low == 1:
                    return low
                if L[mid] <= a:
                    return bSearch(L, a, mid, high)
                else:
                    return bSearch(L, a, low, mid)

        return bSearch(L, a, s, e)

    overlap_idx = []
    seq_linenum_dic = dict()

    for idx in range(0, df1.index.max() + 1):
        seqname = str(df1.loc[idx, "Seq"])
        seq_linenum_dic[seqname] = seq_linenum_dic.get(seqname, []) + [idx]

    for linenum in seq_linenum_dic.values():
        if len(linenum) == 1:
            continue
        else:
            s_lst = []
            e_lst = []
            for lidx in linenum:
                s_lst.append(df1.loc[lidx, "Start"])
                e_lst.append(df1.loc[lidx, "End"])
            for nidx in range(len(linenum) - 1):
                i = Search(s_lst, e_lst[nidx], nidx + 1, len(linenum) - 1)
                if i != None:
                    overlap_idx.extend(linenum[nidx:i + 1])

    overlap_idx = list(sorted(set(overlap_idx)))
    nonoverlap_idx = list(sorted(set(df1.index) - set(overlap_idx)))

    with open(outputfile1, "w") as f:
        for i in overlap_idx:
            f.write("{:>40}\t{:>10}\t{:8.5f}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
                df1.loc[i, "Seq"], \
                df1.loc[i, "OG"], \
                df1.loc[i, "PI"], \
                df1.loc[i, "Start"], \
                df1.loc[i, "End"], \
                df1.loc[i, "Len"], \
                df1.loc[i, "Node"]))
    with open(outputfile2, "w") as f:
        for i in nonoverlap_idx:
            f.write("{:>40}\t{:>10}\t{:8.5f}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(
                df1.loc[i, "Seq"], \
                df1.loc[i, "OG"], \
                df1.loc[i, "PI"], \
                df1.loc[i, "Start"], \
                df1.loc[i, "End"], \
                df1.loc[i, "Len"], \
                df1.loc[i, "Node"]))

overlap_btw_diff_OGs_unmodified()
print("End Step2.6: find overlap & nonoverlap part between each OG: ", time.asctime(time.localtime(time.time())))