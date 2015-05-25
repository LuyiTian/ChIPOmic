#plot pearson of signalvalue and width
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from parse_metadata import get_data_dir, gene_exp_file
from simulation_len import get_len_num
from peak_regression import get_group
from scipy.stats import pearsonr, linregress
import seaborn as sns
#from regression_gene_exp import get_gene_exp_matrix


def plot_width_sigval():
    width_path = os.path.join(get_data_dir(), "tmp", "H3K4me3_width_peak_pearson.csv")
    sigval_path = os.path.join(get_data_dir(), "tmp", "H3K4me3_signalValue_peak_pearson.csv")
    wDF = pd.read_csv(width_path)
    sDF = pd.read_csv(sigval_path)
    X = wDF["Pearson correlation coefficient"].values
    Y = sDF["Pearson correlation coefficient"].values
    plt.plot(X,Y,'x')
    plt.show()


def get_samplewise():
    path = os.path.join(get_data_dir(), "gene_exp", gene_exp_file)
    exp_DF = pd.read_csv(path, sep='\t')
    exp_DF.drop('E000', axis=1, inplace=True)
    EID_list = exp_DF.columns.tolist()[1:]
    #exp_DF = exp_DF.loc[:,["gene_id",eid]]
    path = os.path.join(get_data_dir(), "tmp", "H3K27me3_signalValue_matrix.csv")
    peak_DF = pd.read_csv(path)
    #peak_DF = peak_DF.loc[:,["gene_id",eid]]
    res = []
    X = None
    Y = None
    exp_gene = exp_DF["gene_id"].tolist()
    for g_id in peak_DF["gene_id"].tolist():
        if g_id in exp_gene:
            if X == None:
                X = np.log(peak_DF.loc[peak_DF["gene_id"]==g_id,EID_list].values+1.)
            else:
                X = np.vstack((X,np.log(peak_DF.loc[peak_DF["gene_id"]==g_id,EID_list].values+1.)))
            if Y == None:
                Y = np.log(exp_DF.loc[exp_DF["gene_id"]==g_id,EID_list].values+1.)
            else:
                Y = np.vstack((Y,np.log(exp_DF.loc[exp_DF["gene_id"]==g_id,EID_list].values+1.)))
    for i, eid in enumerate(EID_list):
        print i
        res.append(pearsonr(X[:,i],Y[:,i])[0])
        print pearsonr(X[:,i],Y[:,i])
    #plt.plot(X,Y,'x')
    #print linregress(X,Y)
    dat = pd.DataFrame({'EID':EID_list,'correlation':res})
    dat.to_csv(os.path.join(get_data_dir(), "tmp", "H3K27me3_cor_signalValue_sample.csv"), index=False)


def plot_sample_exp():
    len_dict, _ = get_len_num("H3K4me3")
    DF = pd.read_csv(os.path.join(get_data_dir(), "tmp", "cor_signalValue_sample.csv"))
    X = DF['correlation'].values
    print X
    Y = np.array([np.median(len_dict[eid]) for eid in DF['EID'].tolist()])
    print Y
    dat = pd.DataFrame({'correlation coefficient':X,'median of width':Y})
    sns.lmplot('correlation coefficient','median of width',data=dat)
    plt.show()
    print linregress(X,Y)


def plot_sample_exp_group():
    def get_90q(arr):
        arr.sort()
        return arr[int(0.9*len(arr))]
    len_dict, _ = get_len_num("H3K4me3")
    eid_dic, c_dic = get_group()
    DF = pd.read_csv(os.path.join(get_data_dir(), "tmp", "cor_signalValue_sample.csv"))
    E = DF['EID'].tolist()
    label = []
    width = 0.35
    for ith, key in enumerate(eid_dic):
        X = DF.loc[DF['EID'].isin(eid_dic[key]),'correlation'].values
        Y = np.array([get_90q(len_dict[eid]) for eid in eid_dic[key]])# if eid in E])
        #plt.scatter(X, Y, c=c_dic[key],label=key)
        plt.bar(ith, np.mean(Y), width=width, color=c_dic[key])#, yerr=max(np.std(X),.1))
        label.append(key)
    #plt.legend()
    plt.xticks(np.arange(len(eid_dic))+width/2., label,rotation='vertical')
    plt.show()



if __name__ == '__main__':
    #get_samplewise()
    plot_sample_exp_group()