'''
Created on 2015-5-18

@author: i7
'''
from parse_metadata import get_data_dir,get_full_EID_list
import pandas as pd
from pandas import Series
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot


def find_top(mark,cut_off=0.9):
    a = ['chr'+str(n) for n in range(1,23)]
    a.append('chrX')
    a.append('chrY')
    DF = pd.DataFrame(columns=['chromStart','chromEnd','signalValue','EID','chrom'])
    for chrom in a:
        path = os.path.join(get_data_dir(), "tmp", "{0}-{1}.csv".format(chrom, mark));
        DF_tmp = pd.read_csv(path, sep='\t')
        DF_tmp['chrom'] = chrom
        DF = pd.DataFrame.append(DF,DF_tmp)
    DF['length'] = DF.loc[:, 'chromEnd'].values-DF.loc[:, 'chromStart'].values
    DF = DF.sort('length',ascending=False)
    DF = DF.values[0:len(DF.values)*cut_off]
    DF = pd.DataFrame(DF,columns=['chromStart','chromEnd','signalValue','EID','chrom','length'])
    # DF.to_csv(os.path.join(get_data_dir(), "tmp", "top width -{0}.csv".format(cut_off)), 
    #          sep='\t', index=False)
    print DF
    
def distribution(cut_off = 0.01):
    path = os.path.join(get_data_dir(), "tmp", "top width -{0}.csv".format(cut_off));
    DF = pd.read_csv(path, sep='\t')
    full_EID = get_full_EID_list()
    print DF
    number = []
    for id in full_EID:
        number.append(len(DF.loc[DF['EID'] == id]))
    f, ax = plt.subplots(figsize=(30, 20))
    data = Series(number,index = full_EID)
    data.plot(kind='bar',alpha=0.7)
    plt.show()
    
def percentage():
    path1 = os.path.join(get_data_dir(), "tmp", "top width -0.01.csv");
    DF1 = pd.read_csv(path1, sep='\t')
    path2 = os.path.join(get_data_dir(), "tmp", "top width -1.csv");
    DF2 = pd.read_csv(path2, sep='\t')
    full_EID = get_full_EID_list()
    ratio = [];
    for id in full_EID:
        ratio.append((len(DF1.loc[DF1['EID'] == id])+0.0)/len(DF2.loc[DF2['EID'] == id]))
    f, ax = plt.subplots(figsize=(30, 20))
    data = Series(ratio,index = full_EID)
    data.plot(kind='bar',alpha=0.7)
    plt.show()
    
def all_distribution():
    path = os.path.join(get_data_dir(), "tmp", "top width -1.csv");
    DF = pd.read_csv(path, sep='\t')
    col = ["10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"]
    EID = get_full_EID_list()
    cut_off = [100000,1877,1136,731,482,384,280,233,197,170,1]
    data = {}
    for i in col:
        data[i] = [];
        
    for i in EID:
        DF_tmp = DF.loc[(DF["EID"] == i)]
        for k in range(0,10,1):
            DF_tmp_tmp = DF_tmp.loc[(DF["length"] < cut_off[k]) & (cut_off[k+1] <= DF["length"])]
            data[col[k]].append((len(DF_tmp_tmp)+0.0)/len(DF_tmp))
        
    DF_after = pd.DataFrame(data,columns=col,index=EID)
    print DF_after
    DF_after.to_csv(os.path.join(get_data_dir(),'tmp','all_distribution.csv'),sep='\t')
    

if __name__ == '__main__':
    '''
    mark = 'H3K4me3'
    find_top(mark)
    '''
    # distribution()
    # percentage()
    all_distribution()
    
