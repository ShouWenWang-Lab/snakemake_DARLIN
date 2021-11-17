import argparse
import os

import cospar as cs
import numpy as np
import pandas as pd
from scipy.io import loadmat

import carlinhf as hf

parser = argparse.ArgumentParser(description="Generate heatmap etc")
parser.add_argument(
    "--data_path",
    type=str,
    default=".",
    help="Path of the CARLIN data output",
)

parser.add_argument(
    "--SampleList",
    type=str,
    default="",
    help="A string of sample names, separated by comma",
)

data_path = parser.parse_args().data_path
SampleList= parser.parse_args().SampleList.split(',')

def config(data_path):
    cs.settings.data_path = os.path.join(data_path, "merge_all")
    cs.settings.figure_path = os.path.join(data_path, "merge_all")
    cs.settings.verbosity = 0  # range: 0 (error),1 (warning),2 (info),3 (hint).
    cs.settings.set_figure_params(
        format="png", figsize=[4, 3.5], dpi=300, fontsize=14, pointsize=3, dpi_save=300
    )
    cs.hf.set_up_folders()  # setup the data_path and figure_path


config(data_path)

tmp_list = []
for sample in SampleList:
    base_dir = os.path.join(data_path, f"{sample}")
    df_tmp = hf.load_allele_info(base_dir)
    df_tmp["sample"] = sample.split("_")[0]
    df_tmp["mouse"] = sample.split("-")[0]
    tmp_list.append(df_tmp)
df_all_0 = pd.concat(tmp_list).rename(columns={"UMI_count": "obs_UMI_count"})
adata_orig = hf.generate_adata(df_all_0)

adata_sub = adata_orig[:, adata_orig.uns["multicell_clones"]][
    adata_orig.obs["cells_from_multicell_clone"]
]
adata_sub.obsm["X_clone"] = adata_sub.X
if (adata_sub.X.shape[0]<3) or  (adata_sub.X.shape[1]<2):
    print("No informative clonal data. Skipped")
else:
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
                adata_sub[adata_sub.obs["sample"] == xx].obsm["X_clone"].sum(0).A.flatten()
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

    adata_sub.obs["state_info"] = adata_sub.obs["cell_id"]
    adata_sub.uns["data_des"] = ["refined"]

    cs.pl.barcode_heatmap(
        adata_sub,
        selected_times="0",
        color_bar=True,
        fig_height=10,
        fig_width=10,
        y_ticks=None,  # adata_sub.var_names,
        x_label="Allele",
        y_label="Mutation",
    )

    cs.tl.fate_coupling(adata_sub, selected_times="0", source="X_clone", method="SW")
    cs.pl.fate_coupling(adata_sub, source="X_clone")

    df_info.to_csv(os.path.join(data_path, "merge_all", "X_clone_info.csv"))
