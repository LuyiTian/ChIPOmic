#find correlation between averaged sample peak width and gene expression
import pandas as pd
import numpy as np
import os
from parse_metadata import get_data_dir
from simulation_len import get_len_num
from regression_gene_exp import get_gene_exp_matrix

from scipy.stats import pearsonr

import matplotlib.pyplot as plt
import seaborn as sns


def width_exp_r(mark="H3K4me3"):
    def get_90quantile(arr):
        arr.sort()
        return arr[int(0.9*len(arr))]
    gene_id, EID_list, exp_matrix = get_gene_exp_matrix()
    len_dict, _ = get_len_num(mark)
    quantile_arr = np.array([get_90quantile(len_dict[EID]) for EID in EID_list])
    result = []
    for ith, g_id in enumerate(gene_id):
        tmp = pearsonr(quantile_arr, exp_matrix[ith,:])
        result.append((g_id,tmp[0],tmp[1]))
    result = pd.DataFrame(result,
    columns=['gene_id','Pearson correlation coefficient','2-tailed p-value'])
    result = result.sort("Pearson correlation coefficient")
    for i in result["gene_id"].values[-200:]:
        print i
    result.to_csv(os.path.join(get_data_dir(), "tmp", "{0}_avg_width_exp.csv".format(mark)), index=False)

def plot_avg_width_exp(mark="H3K4me3"):
    def get_90quantile(arr):
        arr.sort()
        return arr[int(0.9*len(arr))]
    gene_id, EID_list, exp_matrix = get_gene_exp_matrix()
    _, len_dict = get_len_num(mark)
    quantile_arr = np.array([get_90quantile(len_dict[EID]) for EID in EID_list])
    gene_avg = np.mean(exp_matrix, axis=0)
    print quantile_arr
    print gene_avg
    result = pd.DataFrame({'quantile90':quantile_arr,"gene_avg":gene_avg})
    sns.lmplot('quantile90','gene_avg',result)
    plt.show()

if __name__ == '__main__':
    plot_avg_width_exp()
