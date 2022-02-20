import skbio.io
import pandas as pd
import concurrent.futures
import time
import sys

print("Start Step2.2: ", time.asctime(time.localtime(time.time())))

inputfile1=sys.argv[1]#"inDir/database/COGFC_GN.tsv.gz"
inputfile2=sys.argv[2]#"runDir/DIAMOND/DIAMOND_raw_output.m8" 
outputfile=sys.argv[3]#"runDir/DIAMOND/qseq_sseq_GNs_num.tsv"

COGFC2GN_df = pd.read_csv(inputfile1, \
                        sep="\t",\
                        names=["COGFC", "GN"],\
                        compression="gzip")

diamond_df = skbio.io.read(inputfile2,\
                        format="blast+6", \
                        into=pd.DataFrame, default_columns=True)

def find_OG_for_sseqid(diamond_row_idx):
    idx = diamond_row_idx
    sseqid = diamond_df.loc[idx, "sseqid"]
    sseqid_GN_df = COGFC2GN_df[(COGFC2GN_df["COGFC"] == str(sseqid))]
    
    GN_list = [i for i in sseqid_GN_df["GN"]]
    GN_num = len(GN_list)
    GN_str = ",".join(GN_list)
    
    new_df = pd.DataFrame({"qseqid": diamond_df.loc[idx, "qseqid"], \
                           "sseqid": diamond_df.loc[idx, "sseqid"], \
                           "GNs": [GN_str], \
                           "GN_num": [GN_num]})
    new_df.to_csv(outputfile, sep="\t", header=False, index=False, mode="a")

with concurrent.futures.ProcessPoolExecutor() as executor:
    row_idx_range = list(range(0, diamond_df.index.max()+1))
    executor.map(find_OG_for_sseqid, row_idx_range)

print("End Step2.2: ", time.asctime(time.localtime(time.time())))