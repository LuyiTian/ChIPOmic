#simulation

import pandas as pd
from parse_metadata import get_full_EID_list, get_data_dir, narrow_peak_col, metadata_filename
import numpy as np
import os
import seaborn as sns


def TSS_len_num(path, mark):
    len_dict = {}
    num_dict = {}
    DF = pd.read_csv(path)
    for EID in get_full_EID_list():
        tmp_DF = DF.loc[DF["EID"]==EID]
        len_dict[EID] = np.log(tmp_DF.loc[:, 'chromEnd'].values-tmp_DF.loc[:, 'chromStart'].values)
        num_dict[EID] = len(tmp_DF)
    return len_dict, num_dict

def DNase_len_num():
    DF = pd.read_csv(os.path.join(get_data_dir(), metadata_filename))
    EID_list = DF.loc[DF["MARK"] == 'DNase','EID'].tolist()
    len_dict = {}
    num_dict = {}
    for EID in EID_list:
        path = os.path.join(get_data_dir(), "hm_data", "{0}-{1}.gz".format(EID, 'DNase'))
        DF = pd.read_csv(path, sep='\t', compression='gzip', header=None, names=narrow_peak_col)
        len_dict[EID] = np.log(DF.loc[:, 'chromEnd'].values-DF.loc[:, 'chromStart'].values)
        num_dict[EID] = len(DF)
    return len_dict, num_dict

def get_len_num(mark):
    len_dict = {}
    num_dict = {}
    for EID in get_full_EID_list():
        path = os.path.join(get_data_dir(), "hm_data", "{0}-{1}.gz".format(EID, mark))
        DF = pd.read_csv(path, sep='\t', compression='gzip', header=None, names=narrow_peak_col)
        len_dict[EID] = np.log(DF.loc[:, 'chromEnd'].values-DF.loc[:, 'chromStart'].values)
        num_dict[EID] = DF.loc[:, 'signalValue'].values
    return len_dict, num_dict


def ran_simu(len_dict, num_dict):
    al = None
    for key in len_dict:
        if al == None:
            al = len_dict[key]
        else:
            al = np.hstack((al, len_dict[key]))
    np.random.shuffle(al)
    res = None
    tmp = 0
    for key in num_dict:
        tmp_arr = al[tmp:tmp+num_dict[key]]
        tmp_arr = np.sort(tmp_arr)
        quantile = [int(i*len(tmp_arr)) for i in [0.01,0.05,0.25,0.50,0.75,0.95,0.99]]
        if res == None:
            res = np.array([tmp_arr[it] for it in quantile])
        else:
            res = np.vstack((res, np.array([tmp_arr[it] for it in quantile])))
        tmp = tmp+num_dict[key]
    return res


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    mark = "H3K27me3"
    len_dict, num_dict = get_len_num(mark)
    result = None
    for i in range(100):
        if result == None:
            result = ran_simu(len_dict, num_dict)
        else:
            result = np.vstack((result, ran_simu(len_dict, num_dict)))
    result = pd.DataFrame(result,
        columns=['0.01', '0.05', '0.25', '0.50', '0.75', '0.95', '0.99'])
    result.to_csv(os.path.join(get_data_dir(), "tmp", "{0}_simulation_result.csv".format(mark)),index=False)
    '''
    data = [len_dict[it] for it in get_full_EID_list()[:4]]
    with sns.color_palette("Set2"):
        for d, label in zip(data, get_full_EID_list()[:4]):
            sns.kdeplot(d, shade=False, label=label)
    plt.show()
    '''
