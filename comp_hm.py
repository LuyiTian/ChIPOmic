#compare same histone marks around different tissues
#and find conserved elements
#contributor: Luyi Tian

import os
from parse_metadata import get_data_dir, narrow_peak_col
import pandas as pd
from itertools import izip


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
    res_df.to_csv(os.path.join(get_data_dir(), "tmp", "{0}-{1}.csv".format(chrom, mark)), sep='\t', index=False)


def find_cluster(mark, chrom, max_diff=100):
    path = os.path.join(get_data_dir(), "tmp", "{0}-{1}.csv".format(chrom, mark))
    DF = pd.read_csv(path, sep='\t')
    DF.convert_objects(convert_numeric=True)
    DF['Cluster'] = 0
    current_mid = None
    index_num = 0
    for the_i, row in DF.iterrows():
        if not current_mid:
            DF[the_i, 'Cluster'] = index_num
            current_mid = [(row['chromEnd']+row['chromStart'])*0.5]
        else:
            if (row['chromEnd']+row['chromStart'])*0.5-sum(current_mid)/len(current_mid) < 150:
                DF[the_i, 'Cluster'] = index_num
                current_mid.append((row['chromEnd']+row['chromStart'])*0.5)
            else:
                index_num += 1
                DF[the_i, 'Cluster'] = index_num
                current_mid = [(row['chromEnd']+row['chromStart'])*0.5]
    DF.to_csv(os.path.join(get_data_dir(), "tmp", "{0}-{1}_clustered.csv".format(chrom, mark)), sep='\t', index=False)


if __name__ == '__main__':
    import numpy as np
    EID_list = ['E002', 'E003', 'E004', 'E005', 'E006', 'E007']
    Full_EID_list = ['E'+str(n).zfill(3) for n in range(1, 130)]
    Full_EID_list.remove('E060')
    Full_EID_list.remove('E064')
    mark = 'H3K4me3'
    '''
    path = os.path.join(get_data_dir(), "tmp", "{0}-{1}.csv".format('chr1', mark))
    DF = pd.read_csv(path, sep='\t')
    DF.convert_objects(convert_numeric=True)
    print DF
    tmp = np.array(DF.loc[1:, 'chromStart']) - np.array(DF.loc[:len(DF)-2, 'chromStart'])
    import pylab as pl
    pl.hist([i for i in tmp if i < 2000], bins=200)
    pl.show()
    '''
    find_cluster(mark, 'chr1')
    #organize_by_chrom(Full_EID_list, mark, 'chr1')
    #com_mark(EID_list, mark)
