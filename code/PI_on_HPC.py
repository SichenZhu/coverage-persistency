import pandas as pd
import time
import copy
import concurrent.futures
import numpy as np
import math
from Bio import SeqIO
import random
import sys

inputfile1=sys.argv[1]#"node_pos_cov_info_PI_og_100.tsv"
inputfile2=sys.argv[2]#"scaffolds.fasta"
outputfile_path_TSV=sys.argv[3]#"${runDir}/TPM/TSV"
outputfile_path_TPM=sys.argv[4]#"${runDir}/TPM/"

positionwise_file = inputfile1
contig_file = inputfile2

# load positionwise data from DIAMOND + Bowtie2 output: name: node_pos_cov_info_PI_og_100.tsv in HPC: /SRS301868/TPM
whole_df = pd.read_csv(positionwise_file, sep="\t", usecols=[0, 1, 2, 4, 5], names=["node", "pos", "cov", "pi", "og"],
                       dtype={"node": 'string', "pos": "int", "cov": "float", "pi": "float", "og": "string"}) #,converters={'og': convert_dtype})
whole_df.fillna(value="n", axis=1, inplace=True)

grouped_df = whole_df.groupby(by=["node"])
group_lst = list(grouped_df.groups.keys())

# load scaffolds.fasta file from SPAdes output:
# /scratch/sz2524/2020summer/BIO/SRS301868/SPAdes/12separate/output/scaffolds.fasta
contig_fasta = list(SeqIO.parse(contig_file, "fasta"))

random.seed(7)
random.shuffle(group_lst)


def calculate_TPM_one_contig(df, nodename, DNAseq):
    # print("Node {} received: ".format(nodename), time.asctime(time.localtime(time.time())))
    # record the number of line where the OG changes(include change between NaN to OG
    idx_og_change_lst = [0]  # index of the df: starts from 0 instead of 1

    for idx in range(0, df.index.max()):
        og0 = str(df.loc[idx, "og"])  # current og
        og1 = str(df.loc[idx + 1, "og"])  # next og
        if og0 == og1:
            continue
        else:
            idx_og_change_lst.extend([idx, idx + 1])  # og switches at [idx, idx+1]
    idx_og_change_lst.append(idx + 1)  # end of the index = len(contig) - 1
    
    # print(idx_og_change_lst)  # list is [0, start1, end1, start2, end2,...,last_coord]
    # example (from local pc /SRS301868/TPM/test1117.tsv): [0, 0, 1, 156, 157, 159, 160, 258, 259, 259]

    idx_og_change_lst_copy = copy.deepcopy(idx_og_change_lst)
    name_v = 0
    for j in range(1, len(idx_og_change_lst_copy), 2):
        # j is the odd index in the idx_og_change_lst which represents the end of each annotation
        start, end = idx_og_change_lst_copy[j - 1], idx_og_change_lst_copy[j]  # 0-based
        gene_len = end - start + 1

        if gene_len < 100:
            df.drop([_ for _ in range(start, end + 1)], inplace=True, errors="ignore")
            idx_og_change_lst.remove(start)
            idx_og_change_lst.remove(end)

        else:
            if str(df.loc[end, "og"]).lower() == "n":
                # rename the non-aligned segment (>100 bp) as a new OG
                nog = str(nodename) + "_" + str(name_v)
                name_v += 1
                for n in range(start, end + 1):
                    df.loc[n, "og"] = nog
                    # calculate TPM
                    df.loc[n, "tpm"] = df.loc[n, "cov"] / gene_len
                    df.loc[n, "nt"] = DNAseq[n]  # add the positionwise nucleotide to df
                    df.loc[n, "node"] = nodename

            else:
                # calculate TPM
                for n in range(start, end + 1):
                    df.loc[n, "tpm"] = df.loc[n, "cov"] / gene_len
                    df.loc[n, "nt"] = DNAseq[n]  # add the positionwise nucleotide to df
                    df.loc[n, "node"] = nodename

    # df.reset_index(drop=True, inplace=True)
    # print("Node {} finish TPM: ".format(nodename), time.asctime(time.localtime(time.time())))
    # print(idx_og_change_lst)  # now the list is [start1, end1, start2, end2...]

    df.to_csv("{}/{}_OG_TPM.tsv".format(outputfile_path_TSV, nodename), sep="\t", header=False, index=False,
              columns=["node", "pos", "og", "pi", "tpm", "nt"])
    # print("Node {} finished: ".format(nodename), time.asctime(time.localtime(time.time())))
    return


def helper_TPM_calc(g, gdf=grouped_df, wholeDNA=contig_fasta):
    node = int(g.split("_")[1])
    ndf = gdf.get_group(g).filter(items=["pos", "cov", "pi", "og"])
    ndf.reset_index(drop=True, inplace=True)
    calculate_TPM_one_contig(ndf, node, wholeDNA[node - 1].seq)
    # wholeDNA index + 1 = NODE number; DNAseq here is a list


def concat_diff_df_total_tpm_calculation(cpuname, g_lst=group_lst):  # serial processing here to avoid error
    # get total_df for the whole sample experiment
    og_ttpm_dic = dict()
    for contig_name in g_lst:
        contig_idx = contig_name.split("_")[1]
        df = pd.read_csv("{}/{}_OG_TPM.tsv".format(outputfile_path_TSV, contig_idx), sep="\t", usecols=[0, 1, 2, 3, 4, 5],
                               names=["node", "pos", "og", "pi", "tpm", "nt"])
        gdf = df.groupby(by=["og"])
        for g in gdf.groups.keys():
            subdf = gdf.get_group(g)
            pi = list(subdf["pi"])[0]
            ttpm = subdf["tpm"].sum()
            og_ttpm_dic[(g, pi)] = og_ttpm_dic.get((g, pi), 0) + ttpm
    # print("{} finish all the contigs: ".format(cpuname), time.asctime(time.localtime(time.time())))

    with open("{}/{}temp".format(outputfile_path_TPM, cpuname), "w") as f1:
        for k,v in og_ttpm_dic.items():
            f1.write("{:>10}\t{:8.5f}\t{:8.5f}\n".format(k[0], k[1], v))
    # print("{} finish dictionary output: ".format(cpuname), time.asctime(time.localtime(time.time())))

    total_df = pd.read_csv("{}/{}temp".format(outputfile_path_TPM, cpuname), sep="\t", usecols=[0, 1, 2], names=["og", "pi", "ttpm"])
    total_gdf = total_df.groupby(by=["og"])
    with open("{}/{}.out".format(outputfile_path_TPM, cpuname), "w") as f2:
        for og_name in total_gdf.groups.keys():
            subdf = total_gdf.get_group(og_name)
            pi = list(subdf["pi"])[0]
            ttpm = subdf["ttpm"].sum()
            f2.write("{:>10}\t{:8.5f}\t{:8.5f}\n".format(og_name, pi, ttpm))

    # print("Done {}: ".format(cpuname), time.asctime(time.localtime(time.time())))

print("Start Step6: TPM calculation for each contig:", time.asctime(time.localtime(time.time())))
# deal with positionwise file: "node_pos_cov_info_PI_og_100.tsv"
with concurrent.futures.ProcessPoolExecutor() as executor:
    executor.map(helper_TPM_calc, group_lst)
print("End Step6: TPM calculation for each contig:", time.asctime(time.localtime(time.time())))

print("Start Step6: TPM calculation for whole sample:", time.asctime(time.localtime(time.time())))
# concatenate all the contig and find total TPM per OG
concat_diff_df_total_tpm_calculation("serial")
print("End Step6: TPM calculation for whole sample:", time.asctime(time.localtime(time.time())))
