# peak-peak relationship
import pandas as pd
import os
from parse_metadata import get_data_dir, get_full_EID_list, sample_group_filename
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA, FastICA
from sklearn.lda import LDA


def get_group():
    '''
    Epigenome ID (EID)
    '''
    DF = pd.read_csv(os.path.join(get_data_dir(), sample_group_filename))
    DF = DF.loc[:,('Epigenome ID (EID)', 'GROUP', "COLOR")]
    tmp = DF.groupby('GROUP').groups
    eid_dic = {}
    c_dic = {}
    for key in tmp:
        eid_dic[key] = [DF.loc[i,'Epigenome ID (EID)'] for i in tmp[key]]
        c_dic[key] = DF.loc[tmp[key][0],'COLOR']
    return eid_dic, c_dic


def get_reg_matrix(input_path, to_csv=True, with_gene_id=True, val="signalValue"):
    DF = pd.read_csv(input_path, sep='\t')
    DF["width"] = DF.loc[:, 'chromEnd']-DF.loc[:, 'chromStart']
    EID_list = get_full_EID_list()
    EID = dict([(eid, i) for i,eid in enumerate(EID_list)])
    gene_list = reduce(lambda r, v: v in r and r or r + [v], list(DF["gene_name"]), [])
    gene = dict([(g, i) for i,g in enumerate(gene_list)])
    res_matrix = np.zeros((len(gene), len(EID_list)))
    tmp_gene = DF.loc[0, "gene_name"]
    for row_indexer in xrange(len(DF)):
        tmp_gene = DF.loc[row_indexer,"gene_name"]
        res_matrix[gene[tmp_gene], EID[DF.loc[row_indexer,"EID"]]] =\
            max(res_matrix[gene[tmp_gene], EID[DF.loc[row_indexer,"EID"]]],float(DF.loc[row_indexer,val]))
        if row_indexer % 1000 == 0:
            print row_indexer
    if to_csv:
        result = pd.DataFrame(res_matrix, columns=EID_list)
        if with_gene_id:
            result['gene_id'] = gene_list
        result.to_csv(os.path.join(get_data_dir(), "tmp", "H3K3me3_{0}_matrix.csv".format(val)), index=False)
    return res_matrix


def l_reg(input_path):
    DF = pd.read_csv(input_path)
    DF.drop('gene_id', axis=1, inplace=True)
    #corr_mat = np.corrcoef(DF.as_matrix())
    f, ax = plt.subplots(figsize=(20, 20))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.corrplot(DF.as_matrix().T, annot=False, sig_stars=False,
             diag_names=False, cmap=cmap, ax=ax)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.savefig(os.path.join(get_data_dir(), "tmp", "H3K27me3_corrplot.png"))


def peak_pca(input_path, eid_dic, c_dic):
    EID_list = get_full_EID_list()
    EID = dict([(eid, i) for i,eid in enumerate(EID_list)])
    DF = pd.read_csv(input_path)
    DF.drop('gene_id', axis=1, inplace=True)
    X = DF.as_matrix().T
    pca = PCA(n_components=2)
    X_r = pca.fit(X).transform(X)
    #ica = FastICA(n_components=2)
    #X_r = ica.fit_transform(X)
    print('explained variance ratio (first two components): %s'
      % str(pca.explained_variance_ratio_))
    plt.figure()
    for gp_name in eid_dic:
        ind = [EID[it] for it in eid_dic[gp_name]]
        plt.scatter(X_r[ind,0],X_r[ind,1],label=gp_name, c=c_dic[gp_name])
    plt.legend()
    plt.title('PCA of H3K4me3 signal value around TSS')
    plt.show()

if __name__ == '__main__':
    path = os.path.join(get_data_dir(), "tmp", "H3K4me3 in TSS-40.csv")
    get_reg_matrix(path)
    #path = os.path.join(get_data_dir(), "tmp", "less40_signalValue_matrix.csv")
    #l_reg(path)
    #path = os.path.join(get_data_dir(), "tmp", "signalValue_matrix_reduced.csv")
    #eid_dic, c_dic = get_group()
    #peak_pca(path, eid_dic, c_dic)