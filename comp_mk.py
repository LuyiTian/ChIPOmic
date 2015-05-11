#compute the coralation matrix between marks in different position
import pandas as pd
from parse_metadata import get_data_dir, get_full_EID_list, gene_annotation_filename, gene_annotation_col
import os
import numpy as np


### for plotting
import seaborn as sns
import matplotlib.pyplot as plt
###


def reg_by_chrom(mark, chrom, cut_off=40):
    Full_EID_list = get_full_EID_list()
    Full_EID_dict = dict([(EID, i) for i, EID in enumerate(Full_EID_list)])
    path = os.path.join(get_data_dir(), "tmp", "{0}-{1}_clustered.csv".format(chrom, mark))
    DF = pd.read_csv(path, sep='\t')
    DF.convert_objects(convert_numeric=True)
    cluster_count = DF.loc[:, 'Cluster'].value_counts()
    cluster_count = cluster_count[cluster_count > cut_off]
    cluster_count = cluster_count.to_dict()
    clu_cot_tuple = cluster_count.items()
    clu_cot_tuple.sort(key=lambda x: x[0])
    res_matrix = []
    for ind, _ in clu_cot_tuple:
        sub_df = DF.loc[DF['Cluster'] == ind, ('EID', 'signalValue')]
        tmp = [0.]*len(Full_EID_list)
        for _, row in sub_df.iterrows():
            tmp[Full_EID_dict[row['EID']]] = row['signalValue']
        res_matrix.append(tmp)
    '''
    f, ax = plt.subplots(figsize=(15, 15))
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    sns.corrplot(np.array(res_matrix).T, annot=False, sig_stars=False,
             diag_names=False, cmap=cmap, ax=ax)
    f.tight_layout()
    plt.show()
    '''
    print np.shape(np.array(res_matrix))
    corr_mat = np.corrcoef(np.array(res_matrix))

    tmp1 = np.array(res_matrix)
    x, y = np.shape(tmp1)
    tmp1 = tmp1.flatten()
    np.random.shuffle(tmp1)
    tmp1 = np.reshape(tmp1, (x, y))
    ran_mat1 = np.corrcoef(tmp1)
    tmp = zip(*res_matrix)
    tmp = [np.random.permutation(col) for col in tmp]
    tmp = np.array(tmp).T
    ran_mat = np.corrcoef(tmp)
    sns.distplot(corr_mat.flatten(), hist=False, kde_kws={"shade": True}, label='true data')
    sns.distplot(ran_mat.flatten(), hist=False, kde_kws={"shade": True}, label='permutation')
    sns.distplot(ran_mat1.flatten(), hist=False, kde_kws={"shade": True}, label='random shuffle')
    #plt.hist([corr_mat.flatten(), ran_mat.flatten()], bins=50, label=['{0} on {1}'.format(mark, chrom), 'random permutation'])
    plt.legend()
    plt.show()


def get_TSS(seq_type="protein_coding"):
    '''
    @param:
    @return:
    '''
    path = os.path.join(get_data_dir(), "gene_exp", gene_annotation_filename)
    DF = pd.read_csv(path, sep='\t', header=None, names=gene_annotation_col)
    DF = DF.loc[DF['type'] == seq_type]
    pos_strand = DF.loc[DF['strand'] == 1, ('ENSG_ID', 'chrom', 'chromStart')]
    pos_strand.rename(columns={'chromStart': 'TSS'}, inplace=True)
    neg_strand = DF.loc[DF['strand'] == -1, ('ENSG_ID', 'chrom', 'chromEnd')]
    neg_strand.rename(columns={'chromEnd': 'TSS'}, inplace=True)
    return pd.concat([pos_strand, neg_strand], ignore_index=True)


def com_exp_mk(mark, chrom):
    '''
    @param:
    @return:
    '''
    pass


if __name__ == '__main__':
    #reg_by_chrom('H3K4me3','chr1')
    print get_TSS()
