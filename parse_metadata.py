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
    "luyi": os.path.normpath("/Users/luyi/data/epi_data"),
    "HZL-PC": os.path.normpath("E:/EDUCATION/BINF90002/assign2/epi_data")
}

MARK = ['H3K4me1', 'H3K4me3', 'H3K27me3', 'H3K36me3', 'H3K9me3']

metadata_filename = "jul2013.roadmapData.qc - Consolidated_EpigenomeIDs_QC.csv"

narrow_peak_url = "http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/"

narrow_peak_col = ["chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"]
#the name following the defination of narrow peak file: http://genome.ucsc.edu/FAQ/FAQformat.html#format12

gene_annotation_filename = "Ensembl_v65.Gencode_v10.ENSG.gene_info"
gene_annotation_col = ["ENSG_ID", "chrom", "chromStart", "chromEnd", "strand", "type", "gene symbol", "full name"]

#sample group file
sample_group_filename = "jul2013.roadmapData.qc - Consolidated_EpigenomeIDs_summary_Table.csv"

#gene expression file
gene_exp_file = "57epigenomes.RPKM.csv"
##########


def get_full_EID_list(numerical=False):
    if numerical:
        Full_EID_list = range(1, 130)
        Full_EID_list.remove(60)
        Full_EID_list.remove(64)
    else:
        Full_EID_list = ['E'+str(n).zfill(3) for n in range(1, 130)]
        Full_EID_list.remove('E060')
        Full_EID_list.remove('E064')
    return Full_EID_list


def get_data_dir():
    '''
    LOCATION/ ->
        ./jul2013.roadmapData.qc - Consolidated_EpigenomeIDs_QC.csv  -> metadata
        hm_data/ -> store histone modification data
        tmp/ -> store temporary data during analysis
        eQTL_data/ -> store eQTL data got from http://www.gtexportal.org/home/
        gene_exp/ -> store gene expression data for epigenomic roadmap project
    '''
    try:
        pc_name = os.environ["COMPUTERNAME"]
    except KeyError:
        pc_name = os.environ['USER']
    root_dir = LOCATION[pc_name]
    sub_dir_list = ["hm_data", "tmp", "eQTL_data", "gene_exp"]
    for sub_dir in sub_dir_list:
        if not os.path.isdir(os.path.join(root_dir, sub_dir)):
            # if sub-dictionary not exist, create it.
            os.mkdir(os.path.join(root_dir, sub_dir))

    return root_dir


def dl_narrow_peak(mark=MARK):
    '''
    @param: mark:selected histone modification marks in list form. e.g. ['H3K27me3']
    '''
    DF = pd.read_csv(os.path.join(get_data_dir(), metadata_filename))
    EID_list = sorted(list(set(list(DF.EID))))
    EID_list = ['E083']  # EID_list[80:]  # the file list is too large, just download part of it one time
    for EID in EID_list:
        for m in mark:
            print ('downloading...    ', EID, m)
            urllib.urlretrieve(
                narrow_peak_url+"{0}-{1}.narrowPeak.gz".format(EID, m),
                os.path.join(get_data_dir(), "hm_data", "{0}-{1}.gz".format(EID, m)))
    print ('download finished~~~~~~~~~~~~~~~')
def dl_DNase_peak():
    DF = pd.read_csv(os.path.join(get_data_dir(), metadata_filename))
    EID_list = DF.loc[DF["MARK"] == 'DNase','EID'].tolist()
    for EID in EID_list:
        print ('downloading...    ', EID, 'DNase')
        urllib.urlretrieve(
            narrow_peak_url+"{}-{}.macs2.narrowPeak.gz".format(EID, 'DNase'),
            os.path.join(get_data_dir(), "hm_data", "{0}-{1}.gz".format(EID, 'DNase')))
    print ('download finished~~~~~~~~~~~~~~~')

if __name__ == '__main__':
    '''
    print (get_data_dir())
    print (os.path.join(get_data_dir(), metadata_filename))
    EID, m = 'E002', 'H3K4me1'
    path = os.path.join(get_data_dir(), "hm_data", "{0}-{1}.gz".format(EID, m))
    print (path)
    DF = pd.read_csv(path, sep='\t', compression='gzip', header=None, names=narrow_peak_col)
    import pylab as pl
    #a = DF['chromEnd']-DF['chromStart']
    b = DF['signalValue']
    pl.hist(b, bins=100)
    pl.show()
    #dl_narrow_peak()
    '''
    dl_DNase_peak()
