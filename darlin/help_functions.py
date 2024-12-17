import os
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import cbook, cm, colors
from matplotlib import pyplot as plt
from matplotlib import rcParams
from scipy.io import loadmat
from tqdm import tqdm


def update_CARLIN_dir(CARLIN_root_folder, template):
    if template == "cCARLIN":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/cCARLIN'"
        )
        # os.system(f"cp {CARLIN_root_folder}/cDARLIN/@CARLIN_def/CARLIN_def_cCARLIN.m {CARLIN_root_folder}/cDARLIN/@CARLIN_def/CARLIN_def.m")
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
    elif template == "Tigre_2022_v2":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Tigre_CARLIN_2022_v2'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Tigre_CARLIN_2022_v2"
    elif template == "Rosa":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Rosa_CARLIN'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Rosa_CARLIN"
    elif template == "Rosa_v2":
        os.system(
            f"rsync -avP '{CARLIN_root_folder}/Custom_CARLIN/' '{CARLIN_root_folder}/Rosa_CARLIN_v2'"
        )
        Actual_CARLIN_dir = f"{CARLIN_root_folder}/Rosa_CARLIN_v2"
    else:
        raise ValueError(
            "The input template should be among {Rosa, Tigre_2022, Tigre, cCARLIN, Rosa_v2, Tigre_2022_v2}"
        )
    return os.path.abspath(Actual_CARLIN_dir)


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


def generate_csv(
    data_path: str,
    SampleList: list,
    no_merge_list: list = None,
    plot=True,
    cfg_type="Bulk",
):
    """
    data_path: should be at the level of samples, e.g., path/to/results_read_cutoff_3

    The no_merge_list should be a list of samples that are from negative control. You want them
    to be included in the csv file, but not when generating the merge_all statistics and relevant files.
    """

    if cfg_type.startswith("Bulk"):
        print("------use Bulk----------")
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
            "valid_5_primer",
            "valid_3_primer",
            "valid_2_seq",
            "valid_read_structure",
            "valid_lines",
            "common_UMIs",
            "called_UMIs_total",
        ]

        annotation = [
            "tot_fastq_N",
            "edit_read_fraction",
            "read_threshold",
            "eventful",
            "called",
            "Mean_read_per_edited_tag",
            "Mean_read_per_unedited_tag",
            "edit_tag_fraction",
            "total_alleles",
            "singleton",
            "effective_allele_N",
            "Diversity_index_all",
            "Diversity_index_edited",
            "CARLIN_potential_by_tag",
            "CARLIN_potential_by_allel",
            "valid_5_primer (read_frac)",
            "valid_3_primer (read_frac)",
            "valid_2_seq (read_frac)",
            "valid_read_structure (read_frac)",
            "valid_lines (read_frac)",
            "common_tags (read_frac)",
            "called_tags_total (read_frac)",
        ]
    elif cfg_type.startswith("sc"):
        print("------use single-cell----------")
        selected_fields = [
            "in_fastq:",
            "eventful_CBs_total:",
            "UMI_chosen:",
            "eventful:",
            "called:",
            "Mean reads per edited CB:",
            "Mean reads per unedited CB:",
            "% CBs edited:",
            "Total (including template):",
            "Singletons (including template):",
            "Effective Alleles:",
            "Diversity Index (normalized by all):",
            "Diversity Index (normalized by edited):",
            "Mean CARLIN potential (by UMI):",
            "Mean CARLIN potential (by allele):",
            "valid_5_primer",
            "valid_3_primer",
            "valid_2_seq",
            "valid_read_structure",
            "valid_lines",
            "common_CBs",
            "called_CBs_total",
        ]

        annotation = [
            "tot_fastq_N",
            "edit_read_fraction",
            "read_threshold",
            "eventful",
            "called",
            "Mean_read_per_edited_tag",
            "Mean_read_per_unedited_tag",
            "edit_tag_fraction",
            "total_alleles",
            "singleton",
            "effective_allele_N",
            "Diversity_index_all",
            "Diversity_index_edited",
            "CARLIN_potential_by_tag",
            "CARLIN_potential_by_allel",
            "valid_5_primer (read_frac)",
            "valid_3_primer (read_frac)",
            "valid_2_seq (read_frac)",
            "valid_read_structure (read_frac)",
            "valid_lines (read_frac)",
            "common_tags (read_frac)",
            "called_tags_total (read_frac)",
        ]
    else:
        raise ValueError("Should start with Bulk or sc")

    df_list = []
    for sample in SampleList:
        ## extract from Results
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

        ## add indel information
        if os.path.exists(f"{data_path}/{sample}/indel_freq_vs_length.mat"):
            pooled_data = loadmat(f"{data_path}/{sample}/indel_freq_vs_length.mat")
            ins_freq = pooled_data["ins_freq"][0][0].flatten()
            del_freq = pooled_data["del_freq"][0][0].flatten()
            size_array = np.arange(len(ins_freq)) + 1
            ave_ins_length = np.sum(ins_freq * size_array) / np.sum(ins_freq)
            ave_del_length = np.sum(del_freq * size_array) / np.sum(del_freq)
            ins_del_ratio_ratio_by_eventful_UMI = np.sum(ins_freq) / np.sum(del_freq)

            my_dict["ave_insert_len"] = [ave_ins_length]
            my_dict["ave_del_len"] = [ave_del_length]
            my_dict["ins_del_ratio_ratio_by_eventful_UMI"] = [
                ins_del_ratio_ratio_by_eventful_UMI
            ]
        else:
            my_dict["ave_insert_len"] = [pd.NA]
            my_dict["ave_del_len"] = [pd.NA]
            my_dict["ins_del_ratio_ratio_by_eventful_UMI"] = [pd.NA]
        my_dict["source"] = sample[:2]
        df_temp = pd.DataFrame(my_dict)
        df_list.append(df_temp.set_index("sample"))

    ## merge all
    df_all_0 = pd.concat(df_list).reset_index()
    if no_merge_list is not None:
        SampleList_new = list(set(SampleList).difference(no_merge_list))
        print(f"Use a new list: {SampleList_new}")
    else:
        SampleList_new = SampleList

    df_temp = compute_merge_all_statistics(
        data_path,
        SampleList_new,
        df_all_0,
        default_key="merge_all",
        plot=plot,
    )
    df_all = pd.concat([df_all_0, df_temp])

    os.makedirs(data_path + "/merge_all", exist_ok=True)
    df_all.to_csv(data_path + "/merge_all/refined_results.csv")


def compute_merge_all_statistics(
    data_path,
    SampleList,
    df_sample_csv_0,
    default_key="merge_all",
    plot=True,
):
    """
    generate a merge_all statistics
    """
    df_sample_csv = df_sample_csv_0[df_sample_csv_0["sample"].isin(SampleList)]
    annotation = list(set(df_sample_csv.columns).difference(["sample"]))
    my_dict = {}
    my_dict["sample"] = default_key
    for j, x in enumerate(annotation):
        my_dict[x] = [np.nan]
    annotation_extensive = [
        "tot_fastq_N",
        "eventful",
        "called",
    ]
    annotation_intensive = [
        "edit_read_fraction",
        "read_threshold",
        "Mean_read_per_edited_tag",
        "Mean_read_per_unedited_tag",
        "edit_tag_fraction",
    ]
    for j, x in enumerate(annotation_extensive):
        my_dict[x] = [df_sample_csv[x].sum()]
    for j, x in enumerate(annotation_intensive):
        my_dict[x] = [
            np.sum(np.array(df_sample_csv["tot_fastq_N"]) * np.array(df_sample_csv[x]))
            / df_sample_csv["tot_fastq_N"].sum()
        ]

    df_count = analyze_allele_frequency_count(data_path, SampleList, plot=plot)
    my_dict["singleton"] = df_count.reset_index()["Frequency"].iloc[0]
    my_dict["total_alleles"] = df_count.reset_index()["Frequency"].sum()
    df_all = pd.read_csv(data_path + "/merge_all/allele_UMI_count.csv")
    my_dict["effective_allele_N"] = effective_allele_number(
        df_all[df_all.allele != "[]"]["UMI_count"]
    )
    my_dict["Diversity_index_all"] = my_dict["effective_allele_N"] / my_dict["called"]
    my_dict["Diversity_index_edited"] = (
        my_dict["effective_allele_N"] / my_dict["eventful"]
    )

    df_temp = pd.DataFrame(my_dict)
    return df_temp


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


def load_allele_frequency_over_multiple_samples(data_path: str, SampleList: list):
    """
    data_path: should be at the level of samples, e.g., path/to/results_read_cutoff_3
    """

    df_list = []
    for sample in SampleList:
        df_temp = load_allele_info(os.path.join(data_path, sample))
        df_list.append(df_temp)
        # null_fraction = (
        #     df_temp[df_temp.allele == "[]"]["UMI_count"].iloc[0]
        #     / df_temp["UMI_count"].sum()
        # )
        # print(f"{sample}: {len(df_temp)}; null fraction {null_fraction:.2f}")
        print(f"{sample}: {len(df_temp)}")
    df_raw = pd.concat(df_list).reset_index()
    df_raw["sample_count"] = 1
    df_new = df_raw.groupby("allele", as_index=False).agg(
        {"UMI_count": "sum", "sample_count": "sum"}
    )
    return df_new


def analyze_allele_frequency_count(data_path: str, SampleList: list, plot=True):
    """
    data_path: should be at the level of samples, e.g., path/to/results_read_cutoff_3
    """

    print(f"sample list in generating allele_UMI_count: {SampleList}")
    os.makedirs(data_path + "/merge_all", exist_ok=True)
    df_new = load_allele_frequency_over_multiple_samples(data_path, SampleList)
    df_new.to_csv(f"{data_path}/merge_all/allele_UMI_count.csv")

    df_count = generate_FrequencyCounts(df_new, save_dir=f"{data_path}/merge_all")

    if plot:
        fig, ax = plt.subplots()
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
        f, axs = plt.subplots(
            1, 2, figsize=(8, 4), gridspec_kw=dict(width_ratios=[4, 4])
        )
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
            "eventful",
            "CARLIN_potential_by_tag",
            "edit_tag_fraction",
            "total_alleles",
        ]
        annotation = [""]
        # title='Read threshold=3'

        selected_variables = [
            "eventful",
            "edit_tag_fraction",
            "called",
            "total_alleles",
            "singleton",
            "effective_allele_N",
            "Diversity_index_all",
            "Diversity_index_edited",
            "CARLIN_potential_by_tag",
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
                    # df_1 = df0.loc["merge_all"]
                    # ax.plot(df_1[x_var], df_1[yy], "r*")

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
    result = {}
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
            # axs[0].yaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
            axs[0].set_ylabel("Count")
            axs[0].set_title(sample)

            axs[1].bar(
                [0, 1],
                [count[0] + count[1], count[2] + count[3]],
                tick_label=["A+T", "C+G"],
            )
            axs[1].set_ylabel("Count")
            # axs[1].yaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
            plt.tight_layout()
            plt.savefig(f"{data_path}/{sample}/insertion_pattern.png")
            result[sample] = count

    f, axs = plt.subplots(1, 2, figsize=(8, 4), gridspec_kw=dict(width_ratios=[4, 3]))
    x_label = ["A", "T", "C", "G"]
    count = [np.sum(np.array(tot_insertion_list) == x) for x in x_label]
    result["All"] = count
    axs[0].bar([0, 1, 2, 3], count, tick_label=x_label)
    # axs[0].yaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
    axs[0].set_ylabel("Count")
    axs[0].set_title("All samples")

    axs[1].bar(
        [0, 1], [count[0] + count[1], count[2] + count[3]], tick_label=["A+T", "C+G"]
    )
    axs[1].set_ylabel("Count")
    # axs[1].yaxis.set_major_formatter(mtick.FormatStrFormatter("%.1e"))
    plt.tight_layout()

    os.makedirs(f"{data_path}/merge_all", exist_ok=True)
    plt.savefig(f"{data_path}/merge_all/all_insertion_pattern.png")
    return pd.DataFrame(result, index=x_label)


def plot_cumulative_insert_del_freq(df_input, save_dir):

    freq = df_input["UMI_count"].to_numpy().astype(int)
    allele_annot = df_input["allele"].to_list()
    del_length_all = []
    ins_length_all = []
    substitute_all = []
    tot_insertion_per_allele = []
    tot_deletion_per_allele = []
    ins_seqs_all = ""
    for i, x in enumerate(allele_annot):
        vector_x = x.split(",")
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
    filehandler = open(f"{save_dir}/insertion_deltion_raw_data.pickle", "wb")
    pickle.dump(all_data, filehandler)

    del_length_all = np.array(del_length_all)
    ins_length_all = np.array(ins_length_all)
    substitute_all = np.array(substitute_all)
    # max_L = np.max([np.max(del_length_all), np.max(ins_length_all), 250])
    max_L = 270

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

    fig, ax = plt.subplots(figsize=(4, 3))
    ax = sns.barplot(x=[1, 2, 3, 4], y=nucleotide_hist)
    # ticks = [1, 2, 3, 4]
    # labels = ["A", "T", "C", "G"]
    plt.xticks(
        ticks=[0, 1, 2, 3], labels=["A", "T", "C", "G"]
    )  # Set text labels and properties.
    ax.set_ylabel("Fraction")
    plt.savefig(f"{save_dir}/insertion_nucleotide_distribution_new.pdf")
    plt.savefig(f"{save_dir}/insertion_nucleotide_distribution_new.png", dpi=100)

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
    g.figure.savefig(f"{save_dir}/cumulative_indel_freq.pdf")
    g.figure.savefig(f"{save_dir}/cumulative_indel_freq.png", dpi=100)

    df.to_csv(f"{save_dir}/cumulative_indel_freq.csv")

    df = pd.DataFrame(
        {"Deletion": tot_deletion_per_allele, "Insertion": tot_insertion_per_allele}
    ).melt()
    fig, ax = plt.subplots(figsize=(4, 3))
    if len(df) > 1:
        ax = sns.violinplot(
            data=df, x="variable", y="value", palette=["#d7301f", "#225ea8"]
        )
    ax.set_ylabel("Cumulative edits per allele (bp)")
    ax.set_xlabel("")
    plt.tight_layout()
    fig.savefig(f"{save_dir}/cumulative_edit_length.pdf")
    fig.savefig(f"{save_dir}/cumulative_edit_length.png", dpi=100)


def effective_allele_number(UMI_counts):
    x = np.array(UMI_counts) / np.sum(UMI_counts)
    entropy = -np.sum(np.log2(x) * x)
    return 2**entropy


def run_sbatch(
    command,
    sbatch_mode="intel-sc3",
    mem="10G",
    cores=2,
    time="01:0:0",
    job_name="sbatch",
):
    os.system("mkdir -p log")
    sbatch_command = f'sbatch -p {sbatch_mode} -c {cores} -t {time} --mem={mem} --job-name {job_name} --output=log/{job_name}-%j.o  --error=log/{job_name}-%j.e --mail-type=TIME_LIMIT_90,FAIL,END --wrap="{command}"'
    print(f"submit job:   {sbatch_command}")
    os.system(sbatch_command)


def merge_fastq(data_path_1, data_path_2, data_path_out, SampleList):
    """
    Merge fastq files from different sequencing run
    ```python
    import yaml
    import os
    root_dir_1='/storage/wangshouwenLab/wangshouwen/DATA/DARLIN/20220717_SC_3A'
    with open(f"{root_dir_1}/config.yaml", "r") as stream:
        file = yaml.safe_load(stream)
        SampleList_1 = file["SampleList"]

    root_dir_2='/storage/wangshouwenLab/wangshouwen/DATA/DARLIN/20220717_SC_3A_2'
    with open(f"{root_dir_2}/config.yaml", "r") as stream:
        file = yaml.safe_load(stream)
        SampleList_2 = file["SampleList"]

    SampleList= list(set(SampleList_1 + SampleList_2))
    SampleList_new=sorted(list(set([x.split('_S')[0] for x in SampleList])))
    ```
    """

    for sample in tqdm(SampleList):
        for source in ["R1", "R2"]:
            file_1 = f"{data_path_1}/{sample}_L001_{source}_001.fastq.gz"
            file_2 = f"{data_path_2}/{sample}_L001_{source}_001.fastq.gz"
            file_3 = f"{data_path_out}/{sample}_L001_{source}_001.fastq.gz"
            if os.path.exists(file_1) & os.path.exists(file_2):
                print(f"merge {sample}")
                os.system(f"cat {file_1} {file_2}   > {file_3}")
            elif os.path.exists(file_1) and (not os.path.exists(file_2)):
                print(f"{sample} not exist in data_path_2. Direct copy")
                os.system(f"cp {file_1} {file_3}")
            elif os.path.exists(file_2) and (not os.path.exists(file_1)):
                print(f"{sample} not exist in data_path_1. Direct copy")
                os.system(f"cp {file_2} {file_3}")
            else:
                print(f"{sample} not exist in either data_path_1 or data_path_2")
