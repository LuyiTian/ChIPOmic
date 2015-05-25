#do box plot

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from simulation_len import get_len_num, TSS_len_num, DNase_len_num
from parse_metadata import get_data_dir
from peak_regression import get_group
import os
def get_simulated_pvalue(p=0.001, quantile = '0.75'):
    path = os.path.join(get_data_dir(), "tmp", "H3K4me3_simulation_result.csv")
    DF = pd.read_csv(path)
    res = DF[quantile].values
    res.sort()
    return res[int(p*len(res))], res[int((1-p)*len(res))]



def peak_box(mark="H3K4me3"):
    len_dict, _ = get_len_num(mark)
    median = [(key,np.median(len_dict[key])) for key in len_dict]
    median.sort(key=lambda x:x[1])
    result = [len_dict[eid] for eid,_ in median]
    plt.boxplot(result, labels=[eid for eid,_ in median])
    lo, hi = get_simulated_pvalue()
    plt.plot([i for i in range(len(median))],[lo for i in range(len(median))])
    plt.xticks([i for i in range(len(median))], [eid for eid,_ in median], rotation='vertical')
    plt.show()


def group_peak_box(mark="H3K4me3"):
    eid_dic, c_dic = get_group()
    #len_dict, _ = TSS_len_num(os.path.join(get_data_dir(), "tmp", "H3K4me3 in TSS-40.csv"),mark)
    len_dict,_ = DNase_len_num()
    label = []
    result = []
    for key in eid_dic:
        for eid in eid_dic[key]:
            if eid in len_dict:
                label.append("{0}--{1}".format(eid,key))
                result.append(len_dict[eid])
    plt.boxplot(result)
    #lo, hi = get_simulated_pvalue()
    #plt.plot([i for i in range(len(label))],[lo for i in range(len(label))])
    plt.xticks([i for i in range(len(label))], label, rotation='vertical')
    plt.show()



if __name__ == '__main__':
    group_peak_box("DNase")
    #print get_simulated_pvalue()
