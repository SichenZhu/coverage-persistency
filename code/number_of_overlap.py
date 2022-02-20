import pandas as pd
import numpy as np
import time
import sys

inputfile=sys.argv[1]#overlap_node_s_e_og_pi.tsv"
outputfile1=sys.argv[2]#"record_overlap_interval_len.tsv"
outputfile2=sys.argv[3]#"record_nonverlap_interval_len.tsv"

print("Start Step2.7: choose PI for overlapping segment: ", time.asctime(time.localtime(time.time())))

### find the number of overlap nucleotides

def number_of_overlap():
    df1 = pd.read_csv(inputfile, sep="\t", \
                      names=["Node", "Start", "End", "OG", "PI"])
    # if only want to use the first 3 columns: usecols = [i for i in range(3)]
    def choose_og_pi(og_1, og_2, pi_1, pi_2, len1, len2):
        if og_1 != og_2:  # different og
            if len1 == len2:
                if pi_1 > pi_2:  # choose the og with higher pi
                    og_ = og_1
                    pi_ = pi_1
                elif pi_1 < pi_2:
                    og_ = og_2
                    pi_ = pi_2
                else:
                    og_ = og_1 + "," + og_2
                    pi_ = pi_1
            elif len1 < len2:
                og_ = og_2
                pi_ = pi_2
            elif len2 < len1:
                og_ = og_1
                pi_ = pi_1
        else:
            og_, pi_ = og_1, pi_1
        return og_, pi_

    def overlap_btw_two(t1, t2, new_value_lst):
        order = sorted([t1, t2], key=lambda x: x[0])
        t1, t2 = order[0], order[1]
        s1, e1, og1, pi1, ol1, l1 = t1[0], t1[1], t1[2], t1[3], t1[4], t1[5]
        s2, e2, og2, pi2, ol2, l2 = t2[0], t2[1], t2[2], t2[3], t2[4], t2[5]
        og, pi = choose_og_pi(og1, og2, pi1, pi2, l1, l2)
        # s1 <= s2 after sorting
        if s1 == s2:
            if e2 < e1:
                new_value_lst += [(s1, e1, og1, pi1, 1, l1)]
                # (e2 + 1, e1, og1, pi1, ol1)
            elif e2 == e1:
                new_value_lst += [(s1, e1, og, pi, 1, l1)] # need to choose OG
            else: # e2 > e1
                new_value_lst += [(s2, e2, og2, pi2, 1, l2)]
                # e1 + 1, e2, og2, pi2, 0)
        elif s2 <= e1:  # s1 must smaller than s2 since it is sorted
            if e2 <= e1:
                new_value_lst += [(s1, e1, og1, pi1, 1, l1)]
                #[(s1, s2 - 1, og1, pi1, ol1), \
                #(s2, e2, og, pi, 1), \
                #(e2 + 1, e1, og1, pi1, ol1)]
            #elif e2 == e1:
                #new_value_lst += 
                #[(s1, s2 - 1, og1, pi1, ol1), \
                #(s2, e2, og, pi, 1)]
            else:
                if l1 < l2:
                    new_value_lst += [(s1, s2 - 1, og1, pi1, ol1, int(s2-s1)), \
                                  (s2, e2, og2, pi2, 1, l2)]
                elif l1 == l2:
                    new_value_lst += [(s1, s2 - 1, og1, pi1, ol1, int(s2-s1)), \
                                  (s2, e1, og, pi, 1, int(e1-s2+1)), \
                                  (e1 + 1, e2, og2, pi2, ol2, int(e2-e1))]
                else: # l1 > l2
                    new_value_lst += [(s1, e1, og1, pi1, 1, l1), \
                                      (e1 + 1, e2, og2, pi2, ol2, int(e2-e1))]
        else:  # s2 > e1: doesn't overlap
            # do not overlap (but the protein itself still exists overlap, later or itself)
            new_value_lst += [t1, t2]  # put the original interval back
        return new_value_lst[:-1], new_value_lst[-1]

    overlap_dic = dict()

    for idx in range(0, df1.index.max() + 1):
        seq = int(df1.loc[idx, "Node"])
        s_2 = int(df1.loc[idx, "Start"])
        e_2 = int(df1.loc[idx, "End"])
        og_2 = str(df1.loc[idx, "OG"])
        pi_2 = float(df1.loc[idx, "PI"])
        length2 = e_2 - s_2 + 1
        if seq not in overlap_dic.keys():  # if seq do not appear before; create one
            overlap_dic[seq] = [(s_2, e_2, og_2, pi_2, 0, length2)]  # if overlap then 1; if not then 0
            # since there must be one line that overlaps with first line, I could first put 0 here and later it would be modified
        else:  # seq already exit
            # iterate through all the intervals in overlap_dic[seq]
            next_t = (s_2, e_2, og_2, pi_2, 0, length2)
            nlst = []
            for i in range(len(overlap_dic[seq])):
                nlst, next_t = overlap_btw_two(overlap_dic[seq][i], next_t, nlst)
                if i == len(overlap_dic[seq]) - 1:
                    nlst.extend([next_t])
            overlap_dic[seq] = sorted(nlst, key=lambda x: x[0])

    with open(outputfile1, "w") as f1:
        total_len = 0
        for k, v in overlap_dic.items():
            len_each_seq = 0
            for pair in v:
                s = pair[0]
                e = pair[1]
                pairlen = int(e - s) + 1
                og = pair[2]
                pi = pair[3]
                overlap_or_not = pair[4]

                if overlap_or_not == 1:  # overlap
                    len_each_seq += pairlen
                    f1.write("{:d}\t{:d}\t{:d}\t{:d}\t{:>10}\t{:8.5f}\n".format(k, s, e, pairlen, og, pi))
            total_len += len_each_seq
        f1.write("{:20}\t{:d}".format("Total overlap length", total_len))

    with open(outputfile2, "w") as f2:
        total_len = 0
        for k, v in overlap_dic.items():
            len_each_seq = 0
            for pair in v:
                s = pair[0]
                e = pair[1]
                pairlen = int(e - s) + 1
                og = pair[2]
                pi = pair[3]
                overlap_or_not = pair[4]

                if overlap_or_not == 0:  # nonoverlap
                    len_each_seq += pairlen
                    f2.write("{:d}\t{:d}\t{:d}\t{:d}\t{:>10}\t{:8.5f}\n".format(k, s, e, pairlen, og, pi))
            total_len += len_each_seq
        f2.write("{:20}\t{:d}".format("Total nonoverlap length", total_len))

number_of_overlap()
print("End Step2.7: choose PI for overlapping segment: ", time.asctime(time.localtime(time.time())))
