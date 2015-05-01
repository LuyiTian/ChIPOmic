#compute the coralation matrix between marks in different position
import pandas as pd
from parse_metadata import get_data_dir, get_full_EID_list
import os
import numpy as np


### for plotting
import seaborn as sns
import matplotlib.pyplot as plt
###


def reg_by_chrom(mark, chrom, cut_off=60):
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
    corr_mat = np.corrcoef(np.array(res_matrix))
    plt.hist(corr_mat.flatten(), bins=50)
    plt.show()


if __name__ == '__main__':
    reg_by_chrom('H3K4me3','chr1')