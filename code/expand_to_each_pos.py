import pandas as pd
import numpy as np
import time
import sys

inputfile1=sys.argv[1]#"/scratch/sz2524/2020summer/BIO/SRS301868/Bowtie2/data_analysis/uniq_node_len.tsv"
inputfile2=sys.argv[2]#"./interval_result.tsv"
outputfile1=sys.argv[3]#"./pos_PI.tsv"
outputfile2=sys.argv[4]#"./cov_no_pi.tsv"
outputfile3=sys.argv[5]#"./pi_no_cov.tsv"

print("Start Step5: Expand DIAMOND & Bowtie2 output to each position: ", time.asctime(time.localtime(time.time())))

def expand_to_each_pos():
    node_len = pd.read_csv(inputfile1, sep="\t", \
                           names=["node", "contiglen"])
    nseop = pd.read_csv(inputfile2, sep="\t", \
                        names=["node", "start", "end", "og", "pi"])

    node_seop_dic = dict()
    for n_idx in range(nseop.index.max()):
        node_seop_dic[nseop.loc[n_idx, "node"]] = node_seop_dic.get(nseop.loc[n_idx, "node"], []) \
                                                  + [(int(nseop.loc[n_idx, "start"]), \
                                                      int(nseop.loc[n_idx, "end"]), \
                                                      str(nseop.loc[n_idx, "og"]), \
                                                      "{:8.5f}".format(nseop.loc[n_idx, "pi"]))]

    cov_no_pi_node_lst = []
    pi_no_cov_node_lst = []

    with open(outputfile1, "w") as f:
        for i in range(node_len.index.max() + 1):
            l0 = node_len.loc[i, "contiglen"]
            pos_info = np.zeros((l0,))
            OG_info = np.full((l0,), None)
            PI_info = np.full((l0,), None)

            try:
                for seop in node_seop_dic[int(node_len.loc[i, "node"])]:
                    pos_info[seop[0] - 1:seop[1]] = np.ones((seop[1] - seop[0] + 1,))
                    OG_info[seop[0] - 1:seop[1]] = np.full((seop[1] - seop[0] + 1,), seop[2])
                    PI_info[seop[0] - 1:seop[1]] = np.full((seop[1] - seop[0] + 1,), seop[3])

                for j in range(len(pos_info)):
                    if OG_info[j] != None:
                        f.write("{:6.0f}\t{:s}\t{:>10}\n".format(pos_info[j], PI_info[j], OG_info[j]))
                    else:
                        f.write("{:6.0f}\t{:d}\t{:>3}\n".format(pos_info[j], -1, np.nan))

                del node_seop_dic[int(node_len.loc[i, "node"])]

            ## if contig is in Bowtie2 but not in DIAMOND output (has coverage but no PI)
            except KeyError:
                cov_no_pi_node_lst.append(node_len.loc[i, "node"])
                for j in range(len(pos_info)):
                    f.write("{:6.0f}\t{:d}\t{:>3}\n".format(pos_info[j], -1, np.nan))
    ## if contig is in DIAMOND but not in Bowtie2 output (has PI but no coverage)
    for k in node_seop_dic.keys():
        pi_no_cov_node_lst.append(k)

    with open(outputfile2, "w") as f1:
        with open(outputfile3, "w") as f2:
            for n1 in cov_no_pi_node_lst:
                f1.write("{:d}\n".format(int(n1)))
            for n2 in pi_no_cov_node_lst:
                f2.write("{:d}\n".format(int(n2)))

expand_to_each_pos()
print("Finish Step5: Expand DIAMOND & Bowtie2 output to each position: ", time.asctime(time.localtime(time.time())))