# __author__ = 'HZL'
#match TFBS within peak region of same sample.

import os
from parse_metadata import get_data_dir, narrow_peak_col
import pandas as pd
import numpy as np
#plotting
import seaborn as sns
import matplotlib.pyplot as plt


def get_TFBS():
    '''
    The TFBS data is downloaded at http://deepbase.sysu.edu.cn/chipbase/chipSeq.php
    '''
    path = os.path.join(get_data_dir(), "hm_data", "chipBase_Human_HeLaS3gamma.csv");
    #build a DF with
    TFBS = pd.read_csv(path, sep=',', header=None)

    #TFBS=TFBS.convert_objects(convert_numeric=True)
    data = {'chrom': TFBS[1], 'chromStart': TFBS[2], 'chromEnd': TFBS[3], 'name': TFBS[0]}
    TFBS = pd.DataFrame(data, columns=['chrom', 'chromStart', 'chromEnd', 'name'])
    TFBS.to_csv(os.path.join(get_data_dir(), "tmp", "HeLaS3gamma_TFBS_sorted.csv"),
                sep='\t', index=False)


def toCsv(filename):
    '''
    reorganize data
    '''
    chroms = ['chr' + str(n) for n in range(1, 23)]
    chroms.append('chrX')
    chroms.append('chrY')
    chrom = []
    chromStart = []
    chromEnd = []
    signalValue=[]
    path = os.path.join(get_data_dir(), "hm_data", filename)
    DF = pd.read_csv(path, sep='\t', compression='gzip', header=None, names=narrow_peak_col)
    DF = DF.loc[:, ('chrom', 'chromStart', 'chromEnd', 'signalValue')]
    DF = DF.sort('chrom')
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
    DF_sorted.to_csv(os.path.join(get_data_dir(), "tmp", filename+"-sorted.csv"),
                sep='\t', index=False)

def map_TFBS(EID, mark):
    chroms = ['chr' + str(n) for n in range(1, 23)]
    chroms.append('chrX')
    chroms.append('chrY')
    path = os.path.join(get_data_dir(), "tmp", 'HeLaS3beta_TFBS_sorted.csv')
    TFBS = pd.read_csv(path, sep='\t')
    chrom = []
    chromStart = []
    chromEnd = []
    signalValue = []

    for chrome in chroms:
        path = os.path.join(get_data_dir(), "tmp", "{0}-{1}.gz".format(EID, mark))
        print path
        DF = pd.read_csv(path, sep='\t', compression='gzip', header=None)
        TFBS_tmp = TFBS.loc[TFBS['chrom'] == chrome]
        #select locations fall in
        for index in TFBS_tmp.index:
            tmp = DF.loc[DF[1].values < TFBS_tmp[1].values]
            tmp = tmp.loc[tmp[2].values > TFBS_tmp[2].values]
            for rows in range(0, len(tmp.index), 1):
                chrom.append(chrome)
                chromStart.apped(tmp.chromStart[rows])
                chromEnd.apped(tmp.chromEnd[rows])
                signalValue.append(tmp.signalValue[rows])

        print 'finished-----------', chrome

    data = {'chrom': chrom, 'chromStart': chromeStart, 'EID': EID, 'signalValue': signalValue}
    DF = pd.DataFrame(data, columns=['chrom', 'chromMiddle', 'EID', 'signalValue'])
    DF.to_csv(os.path.join(get_data_dir(), "tmp", "mapped_{0}_{1} .csv".format(mark, EID)),
              sep='\t', index=False)


if __name__ == "__main__":
    EID_ATF = 'E123'  #K562 cell line
    EID_AP2 = 'E117'  #HeLaS3 cell line
    filename="E117-H3K4me1.gz"
    mark = 'H3K4me3'
    #get_TFBS()
    toCsv(filename)
    #map_TFBS(EID_ATF, mark)
