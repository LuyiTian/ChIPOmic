#compute the coralation matrix between marks in different position
import pandas as pd
from parse_metadata import get_data_dir
import os


def reg_by_chrom(mark,chrom, cut_off = 30):
    path = os.path.join(get_data_dir(), "tmp", "{0}-{1}_clustered.csv".format(chrom, mark))
    DF = pd.read_csv(path, sep='\t')
    DF.convert_objects(convert_numeric=True)
    cluster_count = DF.loc[:, 'Cluster'].value_counts()
    cluster_count = cluster_count[cluster_count > cut_off]
    cluster_count = cluster_count.to_dict()
    for key,val in cluster_count.items():
        print key,val
    print len(cluster_count)



if __name__ == '__main__':
    reg_by_chrom('H3K4me3','chr1')