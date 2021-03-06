#regresssion between gene expression and histone marks around gene TSS
import pandas as pd
import numpy as np
from parse_metadata import get_data_dir, gene_exp_file
import os

from scipy.stats import pearsonr

import matplotlib.pyplot as plt

def get_gene_exp_matrix():
    path = os.path.join(get_data_dir(), "gene_exp", gene_exp_file)
    DF = pd.read_csv(path, sep='\t')
    gene_id = list(DF['gene_id'])
    DF.drop('gene_id', axis=1, inplace=True)
    DF.drop('E000', axis=1, inplace=True)
    EID_list = DF.columns.tolist()
    return gene_id, EID_list, DF.as_matrix()


def gene_peak_regression(val='signalValue', to_csv=True):
    gene_id, EID_list, exp_matrix = get_gene_exp_matrix()
    path = os.path.join(get_data_dir(), "tmp", "H3K27me3_{0}_matrix.csv".format(val))
    DF = pd.read_csv(path)
    peak_gene_id = list(DF['gene_id'])
    result = []
    for ith, g_id in enumerate(gene_id):
        print ith,'gene_id',g_id
        if g_id in peak_gene_id:
            tmp = pearsonr(DF.loc[DF["gene_id"] == g_id, EID_list].values[0,:],exp_matrix[ith,:])
            result.append((g_id,tmp[0],tmp[1]))
    if to_csv:
        result = pd.DataFrame(result,
        columns=['gene_id','Pearson correlation coefficient','2-tailed p-value'])
        result.to_csv(os.path.join(get_data_dir(), "tmp", "H3K27me3_{0}_peak_pearson.csv".format(val)),index=False)
    return result


def plot_r_dist():
    path = os.path.join(get_data_dir(), "tmp", "H3K27me3_signalValue_peak_pearson.csv")
    DF = pd.read_csv(path)
    plt.hist(DF["Pearson correlation coefficient"].values, bins=100)
    plt.show()


def reduce_peak_matrix(r_cutoff=0.5, val='signalValue'):
    path = os.path.join(get_data_dir(), "tmp", "less40_signalValue_peak_pearson.csv")
    DF = pd.read_csv(path)
    gene_list = DF.loc[DF["Pearson correlation coefficient"]>r_cutoff,'gene_id'].tolist()
    for i in gene_list:print i
    path = os.path.join(get_data_dir(), "tmp", "{0}_matrix.csv".format(val))
    DF = pd.read_csv(path)
    DF = DF.loc[DF['gene_id'].isin(gene_list)]
    DF.to_csv(os.path.join(get_data_dir(), "tmp", "{0}_matrix_reduced.csv".format(val)), index=False)



if __name__ == '__main__':
    res = gene_peak_regression()
    plot_r_dist()
    #reduce_peak_matrix()
