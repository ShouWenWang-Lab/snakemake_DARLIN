import argparse
import pickle

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.io import loadmat

from darlin import help_functions as hf

hf.set_rcParams(fontsize=16)
sns.set_style("white")

parser = argparse.ArgumentParser()
parser.add_argument(
    "--input_dir",
    type=str,
    default=".",
)
input_dir = parser.parse_args().input_dir

data = loadmat(f"{input_dir}/allele_annotation.mat")
freq = data["allele_freqs"].flatten()[1:]
allele_annot = data["AlleleAnnotation"].flatten()[1:]
del_length_all = []
ins_length_all = []
substitute_all = []
tot_insertion_per_allele = []
tot_deletion_per_allele = []
ins_seqs_all = ""
for i, x in enumerate(allele_annot):
    vector_x = x[0].split(",")
    temp_del = []
    temp_ins = []
    for y in vector_x:
        if "del" in y:
            temp = y.split("del")[0].split("_")
            del_len = int(temp[1]) - int(temp[0])
            temp_del.append(del_len)
            del_length_all += list(np.repeat(del_len, freq[i]))

        if "ins" in y:
            temp = y.split("ins")[1]
            ins_seqs_all += temp * freq[i]  # consider the allele frequency
            temp_ins.append(len(temp))
            ins_length_all += list(np.repeat(len(temp), freq[i]))

        if ">" in y:
            substitute_all += list(np.repeat(1, freq[i]))

    tot_deletion_per_allele += list(np.repeat(np.sum(temp_del), freq[i]))
    tot_insertion_per_allele += list(np.repeat(np.sum(temp_ins), freq[i]))

all_data = {
    "del_length_all": del_length_all,
    "ins_length_all": ins_length_all,
    "substitute_all": substitute_all,
    "tot_deletion_per_allele": tot_deletion_per_allele,
    "tot_insertion_per_allele": tot_insertion_per_allele,
}
filehandler = open(f"{input_dir}/insertion_deltion_raw_data.pickle", "wb")
pickle.dump(all_data, filehandler)

del_length_all = np.array(del_length_all)
ins_length_all = np.array(ins_length_all)
substitute_all = np.array(substitute_all)
max_L = np.max([np.max(del_length_all), np.max(ins_length_all), 250])

del_length_fraction = np.zeros(max_L)
ins_length_fraction = np.zeros(max_L)
for j, x in enumerate(range(max_L)):
    del_length_fraction[j] = np.sum(del_length_all == x)
    ins_length_fraction[j] = np.sum(ins_length_all == x)


del_length_cum = np.cumsum([0] + list(del_length_fraction / np.sum(freq)))
ins_length_cum = np.cumsum([0] + list(ins_length_fraction / np.sum(freq)))
tot = del_length_cum[-1] + ins_length_cum[-1]
del_length_cum = del_length_cum / tot
ins_length_cum = ins_length_cum / tot


nucleotide_hist = []
for x in ["A", "T", "C", "G"]:
    nucleotide_hist.append(np.sum(np.array(list(ins_seqs_all)) == x))
nucleotide_hist = np.array(nucleotide_hist) / np.sum(nucleotide_hist)

ax = sns.barplot(x=[1, 2, 3, 4], y=nucleotide_hist)
# ticks = [1, 2, 3, 4]
# labels = ["A", "T", "C", "G"]
ax.set_xticks([0, 1, 2, 3], ["A", "T", "C", "G"])  # Set text labels and properties.
ax.set_ylabel("Fraction")
plt.savefig(f"{input_dir}/insertion_nucleotide_distribution_new.pdf")


df = pd.DataFrame(
    {
        "Deletion": del_length_cum,
        "Insertion": ins_length_cum,
        "Size": np.arange(len(ins_length_cum)),
    }
)
df1 = df.melt(id_vars="Size", var_name="Mutation type")
g = sns.relplot(data=df1, x="Size", y="value", hue="Mutation type", kind="line")
# g.figure.autofmt_xdate()
g.ax.set_ylabel("Cumulative fraction of edited cells")
g.ax.set_ylim([0, 1])
g.figure.savefig(f"{input_dir}/cumulative_indel_freq.pdf")

df.to_csv(f"{input_dir}/cumulative_indel_freq.csv")


df = pd.DataFrame(
    {"Deletion": tot_deletion_per_allele, "Insertion": tot_insertion_per_allele}
).melt()
fig, ax = plt.subplots(figsize=(4, 3))
ax = sns.violinplot(data=df, x="variable", y="value", palette=["#d7301f", "#225ea8"])
ax.set_ylabel("Cumulative edits per allele (bp)")
ax.set_xlabel("")
plt.tight_layout()
fig.savefig(f"{input_dir}/cumulative_edit_length.pdf")
