#compare same histone marks around different tissues
#and find conserved elements
#contributor: Luyi Tian

import os
from parse_metadata import get_data_dir, narrow_peak_col, get_full_EID_list
import pandas as pd
from itertools import izip
import numpy as np
# import sklearn for mean shift clustering
from sklearn.cluster import MeanShift, estimate_bandwidth


def com_mark(EID_list, mark):
    '''
    @param:
    @return:
    '''
    res_df = pd.DataFrame()
    for EID in EID_list:
        path = os.path.join(get_data_dir(), "hm_data", "{0}-{1}.gz".format(EID, mark))
        print path
        DF = pd.read_csv(path, sep='\t', compression='gzip', header=None, names=narrow_peak_col)
        DF = DF.loc[:, ('chrom', 'chromStart', 'chromEnd', 'signalValue')]
        if res_df.empty:
            res_df = DF
        else:
            tmp_res = pd.DataFrame(columns=['chrom', 'chromStart', 'chromEnd', 'signalValue'])
            for chrom, chr_st, chr_len, vals in\
                    izip(res_df.loc[:, "chrom"], res_df.loc[:, "chromStart"], res_df.loc[:, "chromEnd"]-res_df.loc[:, "chromStart"], res_df.loc[:, 'signalValue']):
                #print '~~~~~~~', chrom, chr_st
                #TODO: find a more idiomatic method
                tmp = DF.loc[(DF["chromStart"]-chr_st < 0.2*chr_len) & (DF["chromStart"]-chr_st > -0.2*chr_len) & (DF["chrom"] == chrom)]
                if not tmp.empty:
                    tmp.loc[tmp['chromStart'] > chr_st, 'chromStart'] = chr_st
                    tmp.loc[tmp['chromEnd'] < chr_st+chr_len, 'chromEnd'] = chr_st+chr_len
                    tmp.loc[:, 'signalValue'] = str(vals)+','+tmp['signalValue'].astype(str)
                    tmp_res = tmp_res.append(tmp, ignore_index=True)
                    #print tmp_res
            res_df = tmp_res
    res_df.to_csv(os.path.join(get_data_dir(), "tmp", "tmp_res_007.csv"), sep='\t')


def organize_by_chrom(EID_list, mark, chrom):
    '''
    @param:
    @return:
    '''
    res_df = pd.DataFrame(columns=['chromStart', 'chromEnd', 'signalValue', 'EID'])
    for EID in EID_list:
        path = os.path.join(get_data_dir(), "hm_data", "{0}-{1}.gz".format(EID, mark))
        print path
        DF = pd.read_csv(path, sep='\t', compression='gzip', header=None, names=narrow_peak_col)
        DF = DF.loc[DF["chrom"] == chrom, ('chromStart', 'chromEnd', 'signalValue')]
        DF['EID'] = EID
        res_df = res_df.append(DF, ignore_index=True)
    res_df = res_df.sort('chromStart')
    res_df.to_csv(os.path.join(get_data_dir(), "tmp", mark,"{0}-{1}.csv".format(chrom, mark)), sep='\t', index=False)


def find_cluster(mark, chrom, max_diff=150):
    '''
    @param:
    @return:
    '''
    path = os.path.join(get_data_dir(), "tmp", "{0}-{1}.csv".format(chrom, mark))
    DF = pd.read_csv(path, sep='\t')
    DF.convert_objects(convert_numeric=True)
    DF['Cluster'] = 0
    current_mid = None
    index_num = 0
    for the_i, row in DF.iterrows():
        if not current_mid:
            DF.ix[the_i, 'Cluster'] = index_num
            current_mid = [(row['chromEnd']+row['chromStart'])*0.5]
        else:
            if -max_diff < (row['chromEnd']+row['chromStart'])*0.5-sum(current_mid)/len(current_mid) < max_diff:
                DF.ix[the_i, 'Cluster'] = index_num
                current_mid.append((row['chromEnd']+row['chromStart'])*0.5)
            else:
                index_num += 1
                DF.ix[the_i, 'Cluster'] = index_num
                current_mid = [(row['chromEnd']+row['chromStart'])*0.5]
    DF.to_csv(os.path.join(get_data_dir(), "tmp", "{0}-{1}_clustered.csv".format(chrom, mark)), sep='\t', index=False)


def Bland_Altman_plot(mark, chrom, EID=None):
    '''
    @param:
    @return:
    perform Bland-Altman_plot on data:
    http://en.wikipedia.org/wiki/Bland%E2%80%93Altman_plot
    '''
    if EID:
        path = os.path.join(get_data_dir(), "hm_data", "{0}-{1}.gz".format(EID, mark))
        DF = pd.read_csv(path, sep='\t', compression='gzip', header=None, names=narrow_peak_col)
        DF = DF.loc[DF["chrom"] == chrom].sort('chromStart')
        DF = DF.loc[(DF["chromStart"] > 180000) & (DF["chromStart"] < 1500000)]
    else:
        path = os.path.join(get_data_dir(), "tmp", "{0}-{1}.csv".format(chrom, mark))
        DF = pd.read_csv(path, sep='\t')
    S_x = 0.5*(DF.loc[:, 'chromEnd'].values+DF.loc[:, 'chromStart'].values)
    S_y = DF.loc[:, 'chromEnd'].values-DF.loc[:, 'chromStart'].values
    import pylab as pl
    pl.plot(S_x, S_y, '.')
    pl.show()


def BA_meanshift_cluster(mark, chrom):
    '''
    @param:
    @return:
    perform mean shift cluster on 2D data:
        ((chromStart+chromEnd)*0.5, chromEnd-chromStart)
    '''
    path = os.path.join(get_data_dir(), "tmp", "{0}-{1}.csv".format(chrom, mark))
    DF = pd.read_csv(path, sep='\t')
    S_x = 0.5*(DF.loc[:, 'chromEnd'].values+DF.loc[:, 'chromStart'].values)
    S_y = DF.loc[:, 'chromEnd'].values-DF.loc[:, 'chromStart'].values
    X = np.hstack((np.atleast_2d(S_x[5000:6000]).T, np.atleast_2d(S_y[5000:6000]).T))
    print X
    bandwidth = estimate_bandwidth(X, quantile=0.1, n_samples=1000)
    ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
    ms.fit(X)
    labels = ms.labels_
    print list(set(labels))
    import matplotlib.pyplot as plt
    from itertools import cycle
    colors = cycle('bgrcmykbgrcmykbgrcmykbgrcmyk')
    for k, col in zip(range(len(list(set(labels)))), colors):
        my_members = labels == k
        plt.plot(X[my_members, 0], X[my_members, 1], col + '.')
    plt.title('Estimated number of clusters: %d' % len(list(set(labels))))
    plt.show()
if __name__ == '__main__':
    '''
    EID_list = ['E002', 'E003', 'E004', 'E005', 'E006', 'E007']
    Full_EID_list = ['E'+str(n).zfill(3) for n in range(1, 130)]
    Full_EID_list.remove('E060')
    Full_EID_list.remove('E064')
    mark = 'H3K4me3'
    path = os.path.join(get_data_dir(), "tmp", "{0}-{1}_clustered.csv".format('chr1', mark))
    DF = pd.read_csv(path, sep='\t')
    DF.convert_objects(convert_numeric=True)
    TT = DF.loc[:, 'Cluster'].value_counts()
    print TT[TT > 50]
    #tmp = np.array(DF.loc[1:, 'chromStart']) - np.array(DF.loc[:len(DF)-2, 'chromStart'])
    import pylab as pl
    print len([i for i in DF.loc[:, 'Cluster'].value_counts() if i > 30])
    pl.hist([i for i in DF.loc[:, 'Cluster'].value_counts() if i > 30], bins=50)
    pl.show()
    '''
    #Bland_Altman_plot('H3K4me3', 'chr1', 'E002')
    #BA_meanshift_cluster('H3K4me3', 'chr1')
    #find_cluster(mark, 'chr1')
    #organize_by_chrom(Full_EID_list, mark, 'chr1')
    #com_mark(EID_list, mark)
    for chrom in ['chr'+str(i) for i in range(1,23)]+['chrX', 'chrY']:
        organize_by_chrom(get_full_EID_list(), "H3K27me3", chrom)
