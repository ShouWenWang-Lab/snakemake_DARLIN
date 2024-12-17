import numpy as np, matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.ticker as mtick
import os

from matplotlib import rcParams
rcParams["figure.dpi"] = 150

import argparse
parser = argparse.ArgumentParser()
parser.add_argument(
    "--input_dir",
    type=str,
    default=".",
)

parser.add_argument(
    "--output_dir",
    type=str,
    default=".",
)

parser.add_argument(
    "--SampleList",
    type=str,
    default="",
)

input_dir = parser.parse_args().input_dir
output_dir = parser.parse_args().output_dir
SampleList = parser.parse_args().SampleList



### Insertion patterns
allele_annotation={}
insertion_annotation={}
tot_insertion_list=[]
for sample in SampleList.split(','):
    f=f'{input_dir}/{sample}/AlleleAnnotations.txt'
    with open(f, 'r', encoding='utf-8') as infile:
        temp_data=infile.readlines()
        allele_list=[xx.split('\n')[0] for xx in temp_data]
        allele_annotation[sample]=allele_list
        insertion_list=[]
        for xx in allele_list:
            for yy in xx.split(','):
                if 'ins' in yy:
                    insertion_list +=list(yy.split('ins')[1])
                    
        tot_insertion_list += insertion_list
        
        f, axs = plt.subplots(1, 2, figsize=(8, 4),gridspec_kw=dict(width_ratios=[4, 3]))
        x_label=['A','T','C','G']
        count=[np.sum(np.array(insertion_list)==x) for x in x_label]
        axs[0].bar([0,1,2,3],count,tick_label=x_label)
        axs[0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        axs[0].set_ylabel('Count')
        axs[0].set_title(sample)


        axs[1].bar([0,1],[count[0]+count[1],count[2]+count[3]],tick_label=['A+T','C+G'])
        axs[1].set_ylabel('Count')
        axs[1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
        plt.tight_layout()
        plt.savefig(f'{output_dir}/{sample}/insertion_pattern.png')

f, axs = plt.subplots(1, 2, figsize=(8, 4),gridspec_kw=dict(width_ratios=[4, 3]))
x_label=['A','T','C','G']
count=[np.sum(np.array(tot_insertion_list)==x) for x in x_label]
axs[0].bar([0,1,2,3],count,tick_label=x_label)
axs[0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
axs[0].set_ylabel('Count')
axs[0].set_title('All samples')


axs[1].bar([0,1],[count[0]+count[1],count[2]+count[3]],tick_label=['A+T','C+G'])
axs[1].set_ylabel('Count')
axs[1].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
plt.tight_layout()

os.makedirs(f'{output_dir}/merge_all',exist_ok=True)
plt.savefig(f'{output_dir}/merge_all/all_insertion_pattern.png')


### allele breakdown by sample
import os
from scipy.io import loadmat
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import rcParams, cm, colors, cbook
rcParams["figure.dpi"] = 150
import numpy as np

f, axs = plt.subplots(1, 2, figsize=(8, 4),gridspec_kw=dict(width_ratios=[4, 3]))
data=loadmat(f'{output_dir}/merge_all/allele_breakdown_by_sample.mat')
allele_summary=(data['allele_breakdown_by_sample']>0).sum(1)
ratio=np.sum(allele_summary==1)/np.sum(allele_summary==2)
ax=sns.histplot(allele_summary,ax=axs[0])
ax.set_ylabel('Allele acount')
ax.set_xlabel('Occurance in # of mouses')
ax.set_title(f'Occ. ratio (1)/(2)={ratio:.1f}')


unique_ratio=np.sum(allele_summary==1)/np.sum(allele_summary>0)
ax=sns.histplot(allele_summary,ax=axs[1])
plt.yscale('log')
ax.set_ylabel('Allele acount')
ax.set_xlabel('Occurance in # of mouses')
ax.set_title(f'Unique fraction={unique_ratio:.3f}')
plt.tight_layout()
plt.savefig(f'{output_dir}/merge_all/allele_breakdown_by_sample.png')


### visualize the data across samples
df=pd.read_csv(f'{output_dir}/merge_all/refined_results.csv',index_col=0)
df_list=[df]

x_var_list=["UMI_eventful",'CARLIN_potential_by_UMI',"edit_UMI_fraction"]
annotation=['']
#title='Read threshold=3'

selected_variables=['UMI_eventful',
                    "edit_UMI_fraction",
                     'UMI_called',
                     'total_alleles',
                     'singleton',
                     'effective_allele_N',
                     'Diversity_index_all',
                     'Diversity_index_edited',
                     'CARLIN_potential_by_UMI',
                     'CARLIN_potential_by_allel',
                     'ins_del_ratio_ratio_by_eventful_UMI',
                     'ave_insert_len',
                     'ave_del_len']


def start_subplot_figure(n_subplots, n_columns=5, fig_width=14, row_height=3):
    n_rows = int(np.ceil(n_subplots / float(n_columns)))
    fig = plt.figure(figsize=(fig_width, n_rows * row_height))
    return fig, n_rows, n_columns

for x_var in x_var_list:
    fig, n_rows, n_columns=start_subplot_figure(len(selected_variables), 
                            n_columns=4, fig_width=14, row_height=2.5)
    for j, yy in enumerate(selected_variables):
        #fig=plt.figure(figsize=(4,3))
        ax=plt.subplot(n_rows, n_columns,j+1)
        for k,df0 in enumerate(df_list):
            all_samples=list(df0['sample'])
            all_samples.remove('merge_all')
            df0=df0.set_index('sample')
            df_1=df0.loc[all_samples]
            ax.plot(df_1[x_var],df_1[yy],'k^',label=annotation[k])
            df_1=df0.loc['merge_all']
            ax.plot(df_1[x_var],df_1[yy],'r*')
        #ax.legend()
        ax.set_ylabel(yy)
        ax.set_xlabel(x_var)
        #ax.set_title(title)

    plt.tight_layout()
    fig.savefig(f'{output_dir}/merge_all/view_data_across_samples_{x_var}.png')

