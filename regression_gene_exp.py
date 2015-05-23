#regresssion between gene expression and histone marks around gene TSS
import pandas as pd
import numpy as np
from parse_metadata import get_data_dir, gene_exp_file
import os

from scipy.stats import pearsonr

def get_gene_exp_matrix():
    path = os.path.join(get_data_dir(), "gene_exp", gene_exp_file)
    DF = pd.read_csv(path, sep='\t')
    gene_id = list(DF['gene_id'])
    print gene_id[:10]
    DF.drop('gene_id', axis=1, inplace=True)
    DF.drop('E000', axis=1, inplace=True)
    EID_list = DF.columns.tolist()
    return gene_id, EID_list, DF.as_matrix()


def gene_peak_regression(val='signalValue', to_csv=True):
    gene_id, EID_list, exp_matrix = get_gene_exp_matrix()
    path = os.path.join(get_data_dir(), "tmp", "{0}_matrix.csv".format(val))
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
        result.to_csv(os.path.join(get_data_dir(), "tmp", "exp_peak_pearson.csv"),index=False)
    return result


if __name__ == '__main__':
    res = gene_peak_regression()
