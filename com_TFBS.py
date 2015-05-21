# __author__ = 'HZL'
#match TFBS within peak region of same sample.

import os
from parse_metadata import get_data_dir, narrow_peak_col, get_full_EID_list
import pandas as pd
import numpy as np
#plotting
import seaborn as sns
import matplotlib.pyplot as plt


def get_TFBS():
    '''
    The TFBS data is downloaded at http://deepbase.sysu.edu.cn/chipbase/chipSeq.php
    '''
    #read csv file from TFBS file
    path = os.path.join(get_data_dir(), "hm_data", "chipBase_Human_HeLaS3gamma.csv")
    TFBS = pd.read_csv(path, sep=',', header=None)
    #TFBS=TFBS.convert_objects(convert_numeric=True)
    data = {'chrom': TFBS[1], 'chromStart': TFBS[2], 'chromEnd': TFBS[3], 'name': TFBS[0]}
    TFBS = pd.DataFrame(data, columns=['chrom', 'chromStart', 'chromEnd', 'name'])
    #output selected data from dataframe to target file.
    TFBS.to_csv(os.path.join(get_data_dir(), "tmp", "HeLaS3gamma_TFBS_sorted.csv"),
                sep='\t', index=False)


def toCsv(filename):
    '''
    reorganize data by chromStart transition in sorted chrom
    '''
    chroms = ['chr' + str(n) for n in range(1, 23)]
    chroms.append('chrX')
    chroms.append('chrY')
    chrom = []
    chromStart = []
    chromEnd = []
    signalValue=[]
    #load csv file from narrowPeak file
    path = os.path.join(get_data_dir(), "hm_data", filename+".gz")
    DF = pd.read_csv(path, sep='\t', compression='gzip', header=None, names=narrow_peak_col)
    DF = DF.loc[:, ('chrom', 'chromStart', 'chromEnd', 'signalValue')]
    DF = DF.sort('chrom')
    #sort by chromStart in each chromosome
    for chr in chroms:
        DF_tmp=DF.loc[DF["chrom"] == chr].sort('chromStart')
        print 'begin-------',chr
        for j in range(0,len(list(DF_tmp.chrom)),1):
            chrom.append(list(DF_tmp.chrom)[j])
            chromStart.append(list(DF_tmp.chromStart)[j])
            chromEnd.append(list(DF_tmp.chromEnd)[j])
            signalValue.append(list(DF_tmp.signalValue)[j])
        print 'finish-------',chr
    data = {'chrom':chrom,'chromStart':chromStart,'chromEnd':chromEnd,'signalValue':signalValue}
    DF_sorted = pd.DataFrame(data,columns=['chrom','chromStart','chromEnd','signalValue'])
    DF_sorted['width']=DF_sorted['chromEnd'].values-DF_sorted['chromStart']
    #output sorted dataframe to target file.
    DF_sorted.to_csv(os.path.join(get_data_dir(), "tmp", filename+"-sorted.csv"),
                sep='\t', index=False)

def map_TFBS(EID, mark):
    chroms = ['chr' + str(n) for n in range(1, 23)]
    chroms.append('chrX')
    chroms.append('chrY')
    #read csv from TFBS sorted file
    path = os.path.join(get_data_dir(), "tmp", 'HeLaS3gamma_TFBS_sorted.csv')
    TFBS = pd.read_csv(path, sep='\t')
    chrom = []
    chromStart = []
    chromEnd = []
    signalValue = []

    for chrome in chroms:
        #read csv from sorted narrowPeak file
        path = os.path.join(get_data_dir(), "tmp", "{0}-{1}.csv".format(EID, mark))
        DF = pd.read_csv(path, sep='\t', header=None,names=["chrom","chromStart","chromEnd","signalValue"])
        #transform object datatype into long integer datatype.
        DF['chromStart']=DF['chromStart'].astype(long)
        DF['chromEnd']=DF['chromEnd'].astype(long)
        TFBS_tmp = TFBS.loc[TFBS['chrom'] == chrome]
        #select locations fall in peak region
        for start in list(TFBS_tmp['chromStart'].values):
            tmp = DF.loc[DF['chromStart'] < start]
            tmp = tmp.loc[tmp['chromEnd'].values > start+500]
            #add selected values to list
            for rows in range(0, len(tmp.index), 1):
                chrom.append(chrome)
                chromStart.append(list(tmp.chromStart)[rows])
                chromEnd.append(list(tmp.chromEnd)[rows])
                signalValue.append(list(tmp.signalValue)[rows])
        print 'finished-----------', chrome

    data = {'chrom': chrom, 'chromStart': chromStart, 'chromEnd': chromEnd, 'signalValue': signalValue}
    #output
    DF = pd.DataFrame(data, columns=['chrom', 'chromStart', 'chromEnd', 'signalValue'])
    DF.to_csv(os.path.join(get_data_dir(), "tmp", "mapped_{0}_{1} .csv".format(mark, EID)),
              sep='\t', index=False)


if __name__ == "__main__":
    EID_ATF = 'E123'  #K562 cell line
    EID_AP2 = 'E117'  #HeLaS3 cell line
    Full_EID_list=get_full_EID_list()
    mark_list=['H3K4me1', 'H3K4me3'] #['H3K27me3', 'H3K36me3', 'H3K9me3']
    #transform all raw data into organized, sorted csv file
    for EID in Full_EID_list:
        for mark in mark_list:
            filename="{0}-{1}".format(EID,mark)
            toCsv(filename)

    #mark = 'H3K4me1'
    #get_TFBS()
    #map_TFBS(EID_AP2, mark)
