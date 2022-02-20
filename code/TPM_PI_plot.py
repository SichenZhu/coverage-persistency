import sys
import pandas as pd
import numpy as np
import math
import time
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt

inputfile=sys.argv[1] # "serial.out"
outputpdf=sys.argv[2] # "TPM_PI.pdf"

def plot_pi_tpm(filename):
    print("Start Plotting TPM vs. PI figure: ", time.asctime(time.localtime(time.time())))
    df = pd.read_csv("{}".format(filename), sep="\t", usecols=[0, 1, 2], names=["og", "pi", "tpm"])

    def find_corresponding_og2():
        lst = []
        for i in range(df.index.max()):
            if df.loc[i, "pi"] <= 0.1 and df.loc[i, "tpm"] > 20000:  # df.loc[i, "pi"] >= 0.7 and
                # f1.write("{:>10}\t{:8.5f}\t{:8.2f}\n".format(df.loc[i, "og"], df.loc[i, "pi"], df.loc[i, "tpm"]))
                print(df.loc[i])
                lst.append(df.loc[i,"og"].strip())

            # elif df.loc[i, "pi"] <= 0.1 and df.loc[i, "tpm"] > 40000:
                # f2.write("{:>10}\t{:8.5f}\t{:8.2f}\n".format(df.loc[i, "og"], df.loc[i, "pi"], df.loc[i, "tpm"]))
                # print(df.loc[i])
        return lst

    # l = find_corresponding_og2()
    # print(l)

    def plot():
        gb_bins = [-1] + list(np.linspace(0, 1, 11))
        x_label = np.arange(0, 1.2, 0.1)
        labels = list(range(0, 11))
        df["pi_bin"] = pd.cut(df["pi"], bins=gb_bins, labels=labels, include_lowest=True)
        pi_gb = df.groupby(by=["pi_bin"])

        n_lst = []  # n_lst contains all the number in each bin for different pi range
        bins_lst = []
        pi_lst = list(pi_gb.groups.keys())
        tpm_bins_num = 10

        for g in pi_gb.groups.keys():
            try:
                n, bins = np.histogram(list(pi_gb.get_group(g)["tpm"]), bins=tpm_bins_num)  # n, bins for this pi range (i.e. 0.1-0.2)
                n_lst.append(n)
                bins_lst.append(bins)
            except KeyError:
                continue

        plt.figure()
        plt.subplot(111)

        plt.title("TPM vs. PI")
        plt.grid(linestyle="--")
        plt.xlabel("PI")
        plt.ylabel("TPM")
        plt.xticks(ticks=np.arange(0, 1.2, 0.1), labels=["{:3.1f}".format(_) for _ in gb_bins])
        for idx in range(len(n_lst)):
            x = np.full_like(n_lst[idx], fill_value=float(x_label[idx]), dtype=np.double)
            delta = bins_lst[idx][1] - bins_lst[idx][0]
            y = [bins_lst[idx][j] + delta for j in range(tpm_bins_num)]

            for k in range(len(n_lst[idx])):
                if n_lst[idx][k] == 1:
                    plt.scatter(x[k], y[k], s=[math.log(n_lst[idx][k] + 1, 10) * 8], marker="o", alpha=1)
                elif n_lst[idx][k] == 0:
                    continue
                else:
                    plt.scatter(x[k], y[k], s=[math.log(n_lst[idx][k], 10) * 8], marker="o", alpha=1)
        plt.tight_layout()
        plt.savefig(outputpdf, dpi=80)
        print("End Plotting TPM vs. PI figure: ", time.asctime(time.localtime(time.time())))

    plot()

plot_pi_tpm(inputfile)
