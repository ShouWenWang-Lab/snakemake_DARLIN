import os
import sys

import carlinhf as hf
from matplotlib import pyplot as plt

sys.path.append(
    "/Users/shouwenwang/Dropbox (Personal)/shared_folder_LiLi_SW/Code_packages/cospar_master"
)
import cospar_master as cs
import numpy as np
import pandas as pd
from ete3 import AttrFace, NodeStyle, Tree, TreeStyle, faces
from scipy.io import loadmat


def convert_to_tree(parent_map, celltype_names):
    child_map = {
        i: [] for i in set(list(parent_map.values()) + list(parent_map.keys()))
    }
    for i, j in parent_map.items():
        child_map[j].append(i)

    leaf_names = {i: n for i, n in enumerate(celltype_names)}

    def get_newick(n):
        if n in leaf_names:
            return leaf_names[n]
        else:
            return (
                "("
                + ",".join([get_newick(nn) for nn in sorted(child_map[n])[::-1]])
                + ")"
            )

    tree_string = get_newick(np.max(list(child_map.keys()))) + ";"

    t = Tree(tree_string)
    return t


def lineage_analysis_report(
    data_path, SampleList="*", ref_dir=".", data_des="all", color_coding=None
):
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
            format="png",
            figsize=[4, 3.5],
            dpi=300,
            fontsize=14,
            pointsize=3,
            dpi_save=300,
        )
        cs.hf.set_up_folders()  # setup the data_path and figure_path

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
        print(f"Sample: {sample}; allele number: {len(df_tmp)}")
        df_tmp["sample"] = sample.split("_")[0]
        df_tmp["mouse"] = sample.split("-")[0]
        tmp_list.append(df_tmp)
    df_all_0 = pd.concat(tmp_list).rename(columns={"UMI_count": "obs_UMI_count"})
    df_all = hf.query_allele_frequencies(df_ref, df_all_0)
    df_HQ = df_all[df_all.expected_frequency < 10 ** (-4)]
    print("Total unique clone number: {}".format(len(set(df_HQ["allele"]))))
    print("Total cell number: {}".format(len(df_HQ["allele"])))
    adata_orig = hf.generate_adata(df_HQ)

    adata_sub = hf.keep_informative_cell_and_clones(adata_orig)

    if (adata_sub.X.shape[0] < 3) or (adata_sub.X.shape[1] < 2):
        print("No informative clonal data. Skipped")
    else:
        ## generate plots for the coarse-grained data
        # adata_sub.obs["state_info"] = adata_sub.obs["sample"]
        # adata_sub.uns["data_des"] = [f"{data_des}_coarse"]
        # cs.settings.data_path = os.path.join(data_path, "merge_all")
        # cs.settings.figure_path = os.path.join(data_path, "merge_all")
        # cs.pl.barcode_heatmap(
        #     adata_sub, selected_times="0", color_bar=True, fig_height=10, fig_width=16
        # )
        # plt.close()
        # selected_times = "0"
        # cs.settings.set_figure_params(
        #     format="png", figsize=[15, 15], dpi=75, fontsize=8, pointsize=2
        # )
        # coupling, fate_names = cs.pl.fate_coupling_from_clones(
        #     adata_sub, selected_times, method="SW"
        # )
        # plt.close()
        # parent_map, node_mapping = cs.pl.fate_hierarchy_from_clones(
        #     adata_sub, plot_history=False, method="SW"
        # )
        # plt.close()
        # my_tree_coarse = convert_to_tree(parent_map, fate_names)
        # with open(f"{cs.settings.data_path}/{data_des}_coarse_tree.txt", "w") as f:
        #     f.write(my_tree_coarse.write())
        # hf.visualize_tree(
        #     my_tree_coarse,
        #     color_coding=color_coding,
        #     mode="r",
        #     data_des=f"{data_des}_coarse",
        #     figure_path=cs.settings.data_path,
        # )
        # plt.close()

        ## refined heatmap and coupling, no ticks
        adata_sub.obs["state_info"] = adata_sub.obs["cell_id"]
        adata_sub.uns["data_des"] = [f"{data_des}_refined"]
        cs.pl.barcode_heatmap(
            adata_sub, selected_times="0", color_bar=True, fig_height=10, fig_width=16
        )
        plt.close()
        selected_times = "0"
        cs.settings.set_figure_params(
            format="png", figsize=[15, 15], dpi=75, fontsize=8, pointsize=2
        )
        coupling, fate_names = cs.pl.fate_coupling_from_clones(
            adata_sub, selected_times, method="SW"
        )
        plt.close()
        parent_map, node_mapping = cs.pl.fate_hierarchy_from_clones(
            adata_sub, plot_history=False, method="SW"
        )
        plt.close()
        my_tree_refined = convert_to_tree(parent_map, fate_names)
        with open(f"{cs.settings.data_path}/{data_des}_refined_tree.txt", "w") as f:
            f.write(my_tree_refined.write())
        hf.visualize_tree(
            my_tree_refined,
            color_coding=color_coding,
            mode="c",
            data_des=f"{data_des}_refined",
            figure_path=cs.settings.data_path,
            dpi=500,
        )
        plt.close()

        ## get sample information
        cell_N_temp = []
        clone_N_temp = []
        mutation_N_temp = []
        for xx in sorted(set(adata_sub.obs["sample"])):
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
        df_info.to_csv(
            os.path.join(data_path, "merge_all", f"{data_des}_X_clone_info.csv")
        )


if __name__ == "__main__":
    import argparse

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
        default="DATA/CARLIN/20210510_CARLIN_allele_bank/LL_G/CARLIN/results_cutoff_override_3/merge_all",
        help="Dir to allele bank",
    )

    parser.add_argument(
        "--SampleList",
        type=str,
        default="*",
        help="A string of sample names, separated by comma",
    )

    parser.add_argument(
        "--data_des",
        type=str,
        default="all",
        help="A string for file description",
    )

    data_path = parser.parse_args().data_path
    ref_dir = parser.parse_args().ref_dir
    SampleList = parser.parse_args().SampleList
    data_des = parser.parse_args().data_des
    lineage_analysis_report(
        data_path, SampleList=SampleList, ref_dir=ref_dir, data_des=data_des
    )
