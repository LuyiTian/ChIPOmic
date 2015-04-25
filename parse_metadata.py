#parse metadata of epigenomic roadmap data
#and for data retriving, data cleaning and basic data exploration.
#metadata downloaded from here:
#   http://egg2.wustl.edu/roadmap/web_portal/meta.html
#contributor: Luyi Tian
import pandas as pd
import os
import urllib

##########  hard coded parameters
LOCATION = {
    "luyi": os.path.normpath("/Users/luyi/data/epi_data")
}

MARK = ['H3K4me1', 'H3K4me3', 'H3K27me3', 'H3K36me3', 'H3K9me3']

metadata_filename = "jul2013.roadmapData.qc - Consolidated_EpigenomeIDs_QC.csv"

narrow_peak_url = "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/"

narrow_peak_col = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"]
#the name following the defination of narrow peak file: http://genome.ucsc.edu/FAQ/FAQformat.html#format12
##########


def get_data_dir():
    '''
    LOCATION/ ->
        ./jul2013.roadmapData.qc - Consolidated_EpigenomeIDs_QC.csv  -> metadata
        hm_data/ -> store histone modification data
    '''
    try:
        pc_name = os.environ["COMPUTERNAME"]
    except KeyError:
        pc_name = os.environ['USER']
    root_dir = LOCATION[pc_name]
    if not os.path.isdir(os.path.join(root_dir, "hm_data")):
        #if dictionary hm_data/ not exist, create it.
        os.mkdir(os.path.join(root_dir, "hm_data"))
    return root_dir


def dl_narrow_peak(mark=MARK):
    '''
    @param: mark:selected histone modification marks in list form. e.g. ['H3K27me3']
    '''
    DF = pd.read_csv(os.path.join(get_data_dir(), metadata_filename))
    EID_list = sorted(list(set(list(DF.EID))))
    EID_list = EID_list[:10]  # the file list is too large, just download part of it one time
    for EID in EID_list:
        for m in mark:
            print 'downloading...    ', EID, m
            urllib.urlretrieve(
                narrow_peak_url+"{0}-{1}.narrowPeak.gz".format(EID, m),
                os.path.join(get_data_dir(), "hm_data", "{0}-{1}.gz".format(EID, m)))
    print 'download finished~~~~~~~~~~~~~~~'

if __name__ == '__main__':
    print get_data_dir()
    print os.path.join(get_data_dir(), metadata_filename)
    EID, m = 'E002', 'H3K4me1'
    path = os.path.join(get_data_dir(), "hm_data", "{0}-{1}.gz".format(EID, m))
    print path
    DF = pd.read_csv(path, sep='\t', compression='gzip', header=None, names=narrow_peak_col)
    import pylab as pl
    a = list(DF['chromEnd']-DF['chromStart'])
    a = sorted(a)[:50000]
    pl.hist(a, bins=50)
    pl.show()
    #dl_narrow_peak()
