# peak-peak relationship
import pandas as pd
import os
from parse_metadata import get_data_dir, get_full_EID_list
import numpy as np

def get_reg_matrix(input_path, to_csv=True, val="signalValue"):
    DF = pd.read_csv(input_path, sep='\t')
    EID_list = get_full_EID_list()
    EID = dict([(eid, i) for i,eid in enumerate(EID_list)])
    gene= reduce(lambda r, v: v in r and r or r + [v], list(DF["gene_name"]), [])
    gene = dict([(g, i) for i,g in enumerate(gene)])
    res_matrix = np.zeros((len(gene), len(EID_list)))
    tmp_gene = DF.loc[0, "gene_name"]
    for row_indexer in xrange(len(DF)):
        tmp_gene = DF.loc[row_indexer,"gene_name"]
        res_matrix[gene[tmp_gene], EID[DF.loc[row_indexer,"EID"]]] =\
            max(res_matrix[gene[tmp_gene], EID[DF.loc[row_indexer,"EID"]]],float(DF.loc[row_indexer,val]))
        if row_indexer % 1000 == 0:
            print row_indexer
    if to_csv:
        result =  pd.DataFrame(res_matrix, columns=EID_list)
        result.to_csv(os.path.join(get_data_dir(), "tmp", "{0}_matrix.csv".format(val)), index=False)
    return res_matrix


if __name__ == '__main__':
    path = os.path.join(get_data_dir(), "tmp", "H3K4me3_TSS_40.csv")
    get_reg_matrix(path)