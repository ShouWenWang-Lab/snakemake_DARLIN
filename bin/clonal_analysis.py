import argparse
import os

import cospar as cs
import mosaiclineage as hf
import numpy as np
import pandas as pd
from scipy.io import loadmat

parser = argparse.ArgumentParser(description="Generate heatmap etc")
parser.add_argument(
    "--data_path",
    type=str,
    default=".",
    help="Path of the CARLIN data output",
)

parser.add_argument(
    "--ref_dir",
    type=str,
    default=".",
    help="Dir to allele bank",
)

parser.add_argument(
    "--SampleList",
    type=str,
    default="*",
    help="A string of sample names, separated by comma",
)

data_path = parser.parse_args().data_path
ref_dir = parser.parse_args().ref_dir
SampleList = parser.parse_args().SampleList

if SampleList == "*":
    data = loadmat(
        os.path.join(data_path, "merge_all", "allele_breakdown_by_sample.mat")
    )
    SampleList = [xx[0][0] for xx in data["sample_names"]]
else:
    SampleList = SampleList.split(",")


def config(data_path):
    cs.settings.data_path = os.path.join(data_path, "merge_all")
    cs.settings.figure_path = os.path.join(data_path, "merge_all")
    cs.settings.verbosity = 0  # range: 0 (error),1 (warning),2 (info),3 (hint).
    cs.settings.set_figure_params(
        format="png", figsize=[4, 3.5], dpi=300, fontsize=14, pointsize=3, dpi_save=300
    )
    cs.hf.set_up_folders()  # setup the data_path and figure_path


def analysis(adata_temp):
    cs.pl.barcode_heatmap(
        adata_temp,
        color_bar=True,
        fig_height=10,
        fig_width=10,
        y_ticks=None,  # adata_sub.var_names,
        x_label="Allele",
        y_label="Mutation",
        x_ticks=None,
    )
    cs.tl.fate_coupling(adata_temp, selected_times="0", source="X_clone", method="SW")
    cs.pl.fate_coupling(
        adata_temp,
        source="X_clone",
        x_ticks=None,
        y_ticks=None,
        x_label="Allele",
        y_label="Allele",
    )

    cs.tl.fate_hierarchy(adata_temp, source="X_clone")
    my_tree_refined = adata_temp.uns["fate_hierarchy_X_clone"]["tree"]
    with open(f"{cs.settings.data_path}/refined_tree.txt", "w") as f:
        f.write(my_tree_refined.write())
    hf.visualize_tree(
        my_tree_refined,
        color_coding=None,
        mode="c",
        data_des="refined",
        figure_path=cs.settings.data_path,
    )


config(data_path)


df_ref = hf.load_allele_info(ref_dir).rename(columns={"UMI_count": "expect_count"})
tot_count = df_ref["expect_count"].sum()
df_ref["expected_frequency"] = df_ref["expect_count"].apply(lambda x: x / tot_count)
df_ref = df_ref.sort_values("expect_count", ascending=False).filter(
    ["allele", "expected_frequency"]
)

tmp_list = []
for sample in SampleList:
    base_dir = os.path.join(data_path, f"{sample}")
    df_tmp = hf.load_allele_info(base_dir)
    df_tmp["sample"] = sample.split("_")[0]
    df_tmp["mouse"] = sample.split("-")[0]
    tmp_list.append(df_tmp)
df_all_0 = pd.concat(tmp_list).rename(columns={"UMI_count": "obs_UMI_count"})
df_all = hf.query_allele_frequencies(df_ref, df_all_0)
df_HQ = df_all[df_all.expected_frequency < 10 ** (-4)]
print("Clone number: {}".format(len(set(df_HQ["allele"]))))
print("Cell number: {}".format(len(df_HQ["allele"])))
adata_orig = hf.generate_adata(df_HQ)

adata_sub = hf.keep_informative_cell_and_clones(adata_orig)


if (adata_sub.X.shape[0] < 3) or (adata_sub.X.shape[1] < 2):
    print("No informative clonal data. Skipped")
else:
    ## generate plots for the coarse-grained data
    adata_sub.obs["state_info"] = adata_sub.obs["sample"]
    adata_sub.uns["data_des"] = ["coarse"]
    cs.settings.data_path = os.path.join(data_path, "merge_all")
    cs.settings.figure_path = os.path.join(data_path, "merge_all")
    cs.pl.barcode_heatmap(
        adata_sub,
        color_bar=True,
        fig_height=10,
        fig_width=10,
        y_ticks=None,  # adata_sub.var_names,
        x_label="Allele",
        y_label="Mutation",
    )
    cs.tl.fate_coupling(adata_sub, source="X_clone", method="SW")
    cs.pl.fate_coupling(adata_sub, source="X_clone")
    cs.tl.fate_hierarchy(adata_sub, source="X_clone")
    my_tree_coarse = adata_sub.uns["fate_hierarchy_X_clone"]["tree"]
    with open(f"{cs.settings.data_path}/coarse_tree.txt", "w") as f:
        f.write(my_tree_coarse.write())
    hf.visualize_tree(
        my_tree_coarse,
        color_coding=None,
        mode="r",
        data_des="coarse",
        figure_path=cs.settings.data_path,
    )

    adata_sub.obs["state_info"] = adata_sub.obs["cell_id"]
    adata_sub.uns["data_des"] = ["refined"]
    analysis(adata_sub)

    ## get sample information
    cell_N_temp = []
    clone_N_temp = []
    mutation_N_temp = []
    for xx in set(adata_sub.obs["sample"]):
        cell_N_temp.append(np.sum(adata_sub.obs["sample"] == xx))
        clone_N_temp.append(
            len(set(adata_sub[adata_sub.obs["sample"] == xx].obs["allele"]))
        )
        mutation_N_temp.append(
            np.sum(
                adata_sub[adata_sub.obs["sample"] == xx]
                .obsm["X_clone"]
                .sum(0)
                .A.flatten()
                > 0
            )
        )
    df_info = pd.DataFrame(
        {
            "Sample": SampleList,
            "Cell number": cell_N_temp,
            "Clone number": clone_N_temp,
            "mutation number": mutation_N_temp,
        }
    )
    df_info.to_csv(os.path.join(data_path, "merge_all", "X_clone_info.csv"))


for sample in SampleList:
    adata_temp_0 = adata_sub[adata_sub.obs["sample"] == sample.split("_")[0]]
    adata_temp_0.obs["state_info"] = adata_temp_0.obs["cell_id"]
    adata_temp_0.uns["data_des"] = ["refined"]
    adata_temp = hf.keep_informative_cell_and_clones(adata_temp_0)
    cs.settings.data_path = os.path.join(data_path, sample)
    cs.settings.figure_path = os.path.join(data_path, sample)
    analysis(adata_temp)
