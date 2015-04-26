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
    print res_df.to_csv(os.path.join(get_data_dir(), "tmp", "tmp_res_007.csv"), sep='\t')


if __name__ == '__main__':
    EID_list = ['E002', 'E003', 'E004', 'E005', 'E006', 'E007']
    mark = 'H3K4me3'
    com_mark(EID_list, mark)
