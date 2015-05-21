#simulation

import pandas as pd
from parse_metadata import get_full_EID_list, get_data_dir, narrow_peak_col
import numpy as np
import os
import seaborn as sns

def get_len_num(mark):
    len_dict = {}
    num_dict = {}
    for EID in get_full_EID_list():
        path = os.path.join(get_data_dir(), "hm_data", "{0}-{1}.gz".format(EID, mark))
        DF = pd.read_csv(path, sep='\t', compression='gzip', header=None, names=narrow_peak_col)
        len_dict[EID] = np.log(DF.loc[:, 'chromEnd'].values-DF.loc[:, 'chromStart'].values)
        num_dict[EID] = len(DF)
    return len_dict, num_dict


def ran_simu(len_dict,num_dict):
    

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    len_dict, num_dict = get_len_num("H3K4me3")
    al = None
    for key in len_dict:
        if al == None:
            al = len_dict[key]
        else:
            al = np.hstack((al,len_dict[key]))
    data = [len_dict[it] for it in get_full_EID_list()[:4]]+[al]
    with sns.color_palette("Set2"):
        for d, label in zip(data, get_full_EID_list()[:4]+['ALL']):
            sns.kdeplot(d, shade=False, label=label)
    plt.show()
