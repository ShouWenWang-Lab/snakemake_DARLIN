import sys
import pandas as pd

# called_bc_file = "../jobs/20220815_LL837_10X_3A/DARLIN/python_DARLIN/LL837-LF-TA_S6/called_barcodes_by_SW_method.csv"
# carlin_seq_file = "../jobs/20220815_LL837_10X_3A/DARLIN/results_cutoff_override_1/LL837-LF-TA_S6_slim/Actaul_CARLIN_seq.txt"
# allele_annot_file = "../jobs/20220815_LL837_10X_3A/DARLIN/results_cutoff_override_1/LL837-LF-TA_S6_slim/AlleleAnnotations.txt"
# final_file = "../jobs/20220815_LL837_10X_3A/DARLIN/python_DARLIN/LL837-LF-TA_S6/called_barcodes_by_SW_method_allele.csv"
called_bc_file = sys.argv[1]
carlin_seq_file = sys.argv[2]
allele_annot_file = sys.argv[3]
final_file = sys.argv[4]

def hamming_distance(seq1, seq2):
    """
    计算两条相同长度DNA序列的Hamming距离。
    
    参数：
    seq1 (str): 第一条DNA序列。
    seq2 (str): 第二条DNA序列。
    
    返回：
    int: Hamming距离。
    """
    if len(seq1) != len(seq2):
        raise ValueError("DNA序列长度不相等")
    
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

def find_best_match(seq_listA, seq_listB):
    """
    在seq_listA中找到与seq_listB中每个序列最匹配的序列，错配最少。
    
    参数：
    seq_listA (list): 包含较长DNA序列的列表。
    seq_listB (list): 包含较短DNA序列的列表。
    
    返回：
    dict: key为seq_listB中的序列，value为seq_listA中的最佳匹配序列。
    """
    best_matches = {}
    for seqB in seq_listB:
        best_match = None
        min_mismatches = float('inf')
        for seqA in seq_listA:
            if len(seqA) == len(seqB):
                mismatches = hamming_distance(seqA, seqB)
                if mismatches < min_mismatches:
                    min_mismatches = mismatches
                    best_match = seqA
        if best_match is not None:
            best_matches[seqB] = best_match
    return best_matches

# 示例用法
# seq_listA = ["GATTACA", "GACTATA", "GATCACA", "GATTATA"]
# seq_listB = ["GATTACA", "GACTATA"]
# print(find_best_match(seq_listA, seq_listB))  # 输出结果为 {'GATTACA': 'GATTACA', 'GACTATA': 'GACTATA'}

#### Main ####
data_df = pd.read_csv(called_bc_file)
data_df_uniq = data_df[['cell_bc', 'clone_id']].drop_duplicates('clone_id')
clone_barcodes = data_df_uniq['clone_id'].tolist()

## generate new DataFrame with callable information
allele_df = pd.DataFrame({
    "callable_seq": [i.rstrip().replace('-', '') for i in open(carlin_seq_file).readlines()],
    "allele": [i.rstrip() for i in open(allele_annot_file).readlines()]
})
correct_seqs = find_best_match(clone_barcodes, allele_df['callable_seq'])
allele_df['clone_id'] = allele_df['callable_seq'].apply(lambda x: correct_seqs[x])

data_final = pd.merge(data_df, allele_df, how="left", on="clone_id")
data_final.pop('callable_seq')

data_final.to_csv(final_file, index=0)
