import os
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cbook, cm, colors
from matplotlib import pyplot as plt
from matplotlib import rcParams
from scipy.io import loadmat


def update_CARLIN_dir(CARLIN_root_folder, template):
    if template == "cCARLIN":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/cCARLIN'"
        )
        # os.system(f"cp {CARLIN_root_folder}/cCARLIN/@CARLIN_def/CARLIN_def_cCARLIN.m {CARLIN_root_folder}/cCARLIN/@CARLIN_def/CARLIN_def.m")
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/cCARLIN"
    elif template == "Tigre":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Tigre_CARLIN'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Tigre_CARLIN"
    elif template == "Tigre_2022":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Tigre_CARLIN_2022'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Tigre_CARLIN_2022"
    elif template == "Rosa":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Rosa_CARLIN'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Rosa_CARLIN"
    else:
        raise ValueError(
            "The input template should be among {Rosa, Tigre_2022, Tigre, cCARLIN}"
        )
    return Actual_CARLIN_dir


def training_notification(
    to_addr="shouwennotification@gmail.com", msg="pythia training finished"
):
    """
    Note that the msg should not contain ':'.
    """
    import smtplib
    import ssl

    smtp_server = "smtp.gmail.com"
    port = 587  # For starttls
    sender_email = "shouwennotification@gmail.com"
    # to_addr='wangsw09@gmail.com'
    password = "as$3+1245QWcn"  # input("Type your password and press enter: ")

    # msg='pythia training finished---dadadg-=++adgag'

    # Create a secure SSL context
    context = ssl.create_default_context()
    print(msg)

    # Try to log in to server and send email
    try:
        server = smtplib.SMTP(smtp_server, port)
        server.ehlo()  # Can be omitted
        server.starttls(context=context)  # Secure the connection
        server.ehlo()  # Can be omitted
        server.login(sender_email, password)
        server.sendmail(sender_email, to_addr, msg)
        # TODO: Send email here
    except Exception as e:
        # Print any error messages to stdout
        print(e)
    finally:
        server.quit()


def set_rcParams(fontsize=12, color_map=None, frameon=None):
    """Set matplotlib.rcParams to cospar defaults."""
    # check here if you want to customize it: https://matplotlib.org/stable/tutorials/introductory/customizing.html

    # dpi options (mpl default: 100, 100)
    rcParams["figure.dpi"] = 100
    rcParams["savefig.dpi"] = 150

    # figure (mpl default: 0.125, 0.96, 0.15, 0.91)
    rcParams["figure.figsize"] = (6, 4)
    # rcParams["figure.subplot.left"] = 0.18
    # rcParams["figure.subplot.right"] = 0.96
    # rcParams["figure.subplot.bottom"] = 0.15
    # rcParams["figure.subplot.top"] = 0.91

    # lines (defaults:  1.5, 6, 1)
    rcParams["lines.linewidth"] = 1.5  # the line width of the frame
    rcParams["lines.markersize"] = 6
    rcParams["lines.markeredgewidth"] = 1

    # font
    rcParams["font.sans-serif"] = [
        "Arial",
        "Helvetica",
        "DejaVu Sans",
        "Bitstream Vera Sans",
        "sans-serif",
    ]

    fontsize = fontsize
    labelsize = 0.92 * fontsize

    # fonsizes (mpl default: 10, medium, large, medium)
    rcParams["font.size"] = fontsize
    rcParams["legend.fontsize"] = labelsize
    rcParams["axes.titlesize"] = fontsize
    rcParams["axes.labelsize"] = labelsize

    # legend (mpl default: 1, 1, 2, 0.8)
    rcParams["legend.numpoints"] = 1
    rcParams["legend.scatterpoints"] = 1
    rcParams[
        "legend.handlelength"
    ] = 1.5  # change it from 1 to 1.5 to allow seaborn function properly
    rcParams["legend.handletextpad"] = 0.4
    rcParams["pdf.fonttype"] = 42

    # color cycle
    # rcParams["axes.prop_cycle"] = cycler(color=vega_10)

    # axes
    rcParams["axes.linewidth"] = 0.8
    rcParams["axes.edgecolor"] = "black"
    rcParams["axes.facecolor"] = "white"

    # ticks (mpl default: k, k, medium, medium)
    rcParams["xtick.color"] = "k"
    rcParams["ytick.color"] = "k"
    rcParams["xtick.labelsize"] = labelsize
    rcParams["ytick.labelsize"] = labelsize

    # axes grid (mpl default: False, #b0b0b0)
    rcParams["axes.grid"] = False
    rcParams["grid.color"] = ".8"

    # color map
    rcParams["image.cmap"] = "Reds" if color_map is None else color_map

    # spines
    rcParams["axes.spines.right"] = False
    rcParams["axes.spines.top"] = False

    # frame (mpl default: True)
    frameon = False if frameon is None else frameon
    global _frameon
    _frameon = frameon


def generate_csv(data_path: str, SampleList: list):
    """
    data_path: should be at the level of samples, e.g., path/to/results_read_cutoff_3
    """

    selected_fields = [
        "in_fastq:",
        "eventful_UMIs_total:",
        "UMI_chosen:",
        "eventful:",
        "called:",
        "Mean reads per edited UMI:",
        "Mean reads per unedited UMI:",
        "% UMIs edited:",
        "Total (including template):",
        "Singletons (including template):",
        "Effective Alleles:",
        "Diversity Index (normalized by all):",
        "Diversity Index (normalized by edited):",
        "Mean CARLIN potential (by UMI):",
        "Mean CARLIN potential (by allele):",
    ]

    annotation = [
        "tot_fastq_N",
        "edit_read_fraction",
        "read_threshold",
        "UMI_eventful",
        "UMI_called",
        "Mean_read_per_edited_UMI",
        "Mean_read_per_unedited_UMI",
        "edit_UMI_fraction",
        "total_alleles",
        "singleton",
        "effective_allele_N",
        "Diversity_index_all",
        "Diversity_index_edited",
        "CARLIN_potential_by_UMI",
        "CARLIN_potential_by_allel",
    ]
    df_list = []
    for sample in SampleList:
        pooled_data = loadmat(f"{data_path}/{sample}/indel_freq_vs_length.mat")
        ins_freq = pooled_data["ins_freq"][0][0].flatten()
        del_freq = pooled_data["del_freq"][0][0].flatten()
        size_array = np.arange(len(ins_freq)) + 1
        ave_ins_length = np.sum(ins_freq * size_array) / np.sum(ins_freq)
        ave_del_length = np.sum(del_freq * size_array) / np.sum(del_freq)
        ins_del_ratio_ratio_by_eventful_UMI = np.sum(ins_freq) / np.sum(del_freq)

        filename = f"{data_path}/{sample}/Results.txt"
        with open(filename) as file:
            lines = [line.rstrip() for line in file]

        value = []
        for y in selected_fields:
            for x in lines:
                if y in x:
                    temp = [x2 for x2 in x.split(" ") if x2 != ""]
                    if y == "in_fastq:":
                        value.append(temp[-2])
                    else:
                        value.append(temp[-1])

        value = np.array(value).astype(float)
        my_dict = {}
        my_dict["sample"] = sample
        for j, x in enumerate(annotation):
            my_dict[x] = [value[j]]

        my_dict["ave_insert_len"] = [ave_ins_length]
        my_dict["ave_del_len"] = [ave_del_length]
        my_dict["ins_del_ratio_ratio_by_eventful_UMI"] = [
            ins_del_ratio_ratio_by_eventful_UMI
        ]
        my_dict["source"] = sample[:2]
        df_temp = pd.DataFrame(my_dict)
        df_list.append(df_temp.set_index("sample"))

    ## merge all
    df_all_0 = pd.concat(df_list).reset_index()
    my_dict = {}
    my_dict["sample"] = "merge_all"
    for j, x in enumerate(annotation):
        my_dict[x] = [np.nan]
    annotation_extensive = [
        "tot_fastq_N",
        "UMI_eventful",
        "UMI_called",
    ]
    annotation_intensive = [
        "edit_read_fraction",
        "read_threshold",
        "Mean_read_per_edited_UMI",
        "Mean_read_per_unedited_UMI",
        "edit_UMI_fraction",
    ]
    for j, x in enumerate(annotation_extensive):
        my_dict[x] = [df_all_0[x].sum()]
    for j, x in enumerate(annotation_intensive):
        my_dict[x] = [
            np.sum(np.array(df_all_0["tot_fastq_N"]) * np.array(df_all_0[x]))
            / df_all_0["tot_fastq_N"].sum()
        ]

    df_count = analyze_allele_frequency_count(data_path, SampleList)
    my_dict["singleton"] = df_count.reset_index()["Frequency"].iloc[0]
    my_dict["total_alleles"] = df_count.reset_index()["Frequency"].sum()

    df_temp = pd.DataFrame(my_dict)
    df_all = pd.concat([df_all_0, df_temp])

    os.makedirs(data_path + "/merge_all", exist_ok=True)
    df_all.to_csv(data_path + "/merge_all/refined_results.csv")


def load_allele_info(data_path):
    pooled_data = loadmat(os.path.join(data_path, "allele_annotation.mat"))
    allele_freqs = pooled_data["allele_freqs"].flatten()
    alleles = [xx[0][0] for xx in pooled_data["AlleleAnnotation"]]
    return pd.DataFrame({"allele": alleles, "UMI_count": allele_freqs})


def generate_FrequencyCounts(df_raw, save_dir=None):
    """
    df_raw is a pandas object, with
    'allele','UMI_count'

    A speeded up version
    """
    df_input = df_raw.reset_index()
    df_new = df_input.groupby("allele", as_index=False).agg({"UMI_count": "sum"})

    UMI_count = list(df_new["UMI_count"])
    unique_count = np.sort(list(set(UMI_count))).astype(int)
    count_frequency = np.zeros(len(unique_count), dtype=int)
    for j, x in enumerate(unique_count):
        count_frequency[j] = np.sum(UMI_count == x)

    df_count = pd.DataFrame(
        {"UMI_count": unique_count, "Frequency": count_frequency}
    ).set_index("UMI_count")
    if save_dir is not None:
        df_count.to_csv(f"{save_dir}/FrequencyCounts.csv", header=None)
    return df_count


def analyze_allele_frequency_count(data_path: str, SampleList: list):
    """
    data_path: should be at the level of samples, e.g., path/to/results_read_cutoff_3
    """

    os.makedirs(data_path + "/merge_all", exist_ok=True)
    df_list = []
    for sample in SampleList:
        df_temp = load_allele_info(os.path.join(data_path, sample))
        df_list.append(df_temp)
        print(f"{sample}: {len(df_temp)}")
    df_raw = pd.concat(df_list).reset_index()
    df_raw["sample_count"] = 1
    df_new = df_raw.groupby("allele", as_index=False).agg(
        {"UMI_count": "sum", "sample_count": "sum"}
    )
    df_new.to_csv(f"{data_path}/merge_all/allele_UMI_count.csv")

    df_count = generate_FrequencyCounts(df_new, save_dir=f"{data_path}/merge_all")
    ax = sns.scatterplot(
        data=df_count.reset_index(),
        x="UMI_count",
        y="Frequency",
        edgecolor="k",
        color="#d7191c",
    )
    plt.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
    plt.xscale("log")
    plt.xscale("log")
    # ax.set_xlim([0,10.5])

    singleton = df_count.reset_index()["Frequency"].iloc[0]
    singleton_fraction = (
        df_count.reset_index()["Frequency"].iloc[0]
        / df_count.reset_index()["Frequency"].sum()
    )
    ax.set_xlabel("UMI count per allele")
    ax.set_ylabel("Allele number per UMI count")
    ax.set_title(f"Singleton fraction: {singleton_fraction:.2f}")
    plt.tight_layout()
    plt.savefig(f"{data_path}/merge_all/Frequency_count.pdf")

    ## allele breakdown by samples
    f, axs = plt.subplots(1, 2, figsize=(8, 4), gridspec_kw=dict(width_ratios=[4, 4]))
    allele_summary = np.array(df_new["sample_count"])
    ratio = np.sum(allele_summary == 1) / len(df_new)
    ax = sns.histplot(allele_summary, ax=axs[0])
    ax.set_ylabel("Allele acount")
    ax.set_xlabel("Occurance in # of samples")
    ax.set_title(f"Frac. in 1 sample={ratio:.2f}")

    ratio_1_2 = np.sum(allele_summary == 1) / np.sum(allele_summary == 2)
    ax = sns.histplot(allele_summary, ax=axs[1])
    plt.yscale("log")
    ax.set_ylabel("Allele acount")
    ax.set_xlabel("Occurance in # of samples")
    ax.set_title(f"Occu. ratio (1/2)={ratio_1_2:.2f}")
    plt.tight_layout()
    plt.savefig(f"{data_path}/merge_all/allele_breakdown_by_sample.png")
    return df_count


def plot_data_statistics_across_samples(data_path):
    file_path = f"{data_path}/merge_all/refined_results.csv"
    if not os.path.exists(file_path):
        raise ValueError("CSV file not computed yet")
    else:
        df = pd.read_csv(file_path, index_col=0)
        df_list = [df]

        x_var_list = [
            "UMI_eventful",
            "CARLIN_potential_by_UMI",
            "edit_UMI_fraction",
            "total_alleles",
        ]
        annotation = [""]
        # title='Read threshold=3'

        selected_variables = [
            "UMI_eventful",
            "edit_UMI_fraction",
            "UMI_called",
            "total_alleles",
            "singleton",
            "effective_allele_N",
            "Diversity_index_all",
            "Diversity_index_edited",
            "CARLIN_potential_by_UMI",
            "CARLIN_potential_by_allel",
            "ins_del_ratio_ratio_by_eventful_UMI",
            "ave_insert_len",
            "ave_del_len",
        ]

        def start_subplot_figure(n_subplots, n_columns=5, fig_width=14, row_height=3):
            n_rows = int(np.ceil(n_subplots / float(n_columns)))
            fig = plt.figure(figsize=(fig_width, n_rows * row_height))
            return fig, n_rows, n_columns

        for x_var in x_var_list:
            fig, n_rows, n_columns = start_subplot_figure(
                len(selected_variables), n_columns=4, fig_width=14, row_height=2.5
            )
            for j, yy in enumerate(selected_variables):
                # fig=plt.figure(figsize=(4,3))
                ax = plt.subplot(n_rows, n_columns, j + 1)
                for k, df0 in enumerate(df_list):
                    all_samples = list(df0["sample"])
                    all_samples.remove("merge_all")
                    df0 = df0.set_index("sample")
                    df_1 = df0.loc[all_samples]
                    ax.plot(df_1[x_var], df_1[yy], "k^", label=annotation[k])
                    #df_1 = df0.loc["merge_all"]
                    #ax.plot(df_1[x_var], df_1[yy], "r*")

                # ax.legend()
                ax.set_ylabel(yy)
                ax.set_xlabel(x_var)
                # ax.set_title(title)

            plt.tight_layout()
            fig.savefig(f"{data_path}/merge_all/view_data_across_samples_{x_var}.png")


def plot_insertion_patterns(data_path: str, SampleList: list):
    import matplotlib.ticker as mtick

    allele_annotation = {}
    insertion_annotation = {}
    tot_insertion_list = []
    for sample in SampleList:
        f = f"{data_path}/{sample}/AlleleAnnotations.txt"
        with open(f, "r", encoding="utf-8") as infile:
            temp_data = infile.readlines()
            allele_list = [xx.split("\n")[0] for xx in temp_data]
            allele_annotation[sample] = allele_list
            insertion_list = []
            for xx in allele_list:
                for yy in xx.split(","):
                    if "ins" in yy:
                        insertion_list += list(yy.split("ins")[1])

            tot_insertion_list += insertion_list

            f, axs = plt.subplots(
                1, 2, figsize=(8, 4), gridspec_kw=dict(width_ratios=[4, 3])
            )
            x_label = ["A", "T", "C", "G"]
            count = [np.sum(np.array(insertion_list) == x) for x in x_label]
            axs[0].bar([0, 1, 2, 3], count, tick_label=x_label)
            axs[0].yaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
            axs[0].set_ylabel("Count")
            axs[0].set_title(sample)

            axs[1].bar(
                [0, 1],
                [count[0] + count[1], count[2] + count[3]],
                tick_label=["A+T", "C+G"],
            )
            axs[1].set_ylabel("Count")
            axs[1].yaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
            plt.tight_layout()
            plt.savefig(f"{data_path}/{sample}/insertion_pattern.png")

    f, axs = plt.subplots(1, 2, figsize=(8, 4), gridspec_kw=dict(width_ratios=[4, 3]))
    x_label = ["A", "T", "C", "G"]
    count = [np.sum(np.array(tot_insertion_list) == x) for x in x_label]
    axs[0].bar([0, 1, 2, 3], count, tick_label=x_label)
    axs[0].yaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
    axs[0].set_ylabel("Count")
    axs[0].set_title("All samples")

    axs[1].bar(
        [0, 1], [count[0] + count[1], count[2] + count[3]], tick_label=["A+T", "C+G"]
    )
    axs[1].set_ylabel("Count")
    axs[1].yaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
    plt.tight_layout()

    os.makedirs(f"{data_path}/merge_all", exist_ok=True)
    plt.savefig(f"{data_path}/merge_all/all_insertion_pattern.png")
