'''
Created on 2015-5-16

@author: i7
'''
import os
from parse_metadata import get_data_dir,get_full_EID_list
from Bio import SeqIO
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import plot
import math

def paixu(DF):
    a = ['chr'+str(n) for n in range(1,23)]
    a.append('chrX')
    a.append('chrY')
    chrom = [];
    chromStart = [];
    chromEnd = [];
    for i in a:
        DF_tmp = DF.loc[DF["chrom"] == i].sort('chromStart')
        print 'begin',i
        for j in range(0,len(list(DF_tmp.chrom)),1):
            chrom.append(list(DF_tmp.chrom)[j])
            chromStart.append(list(DF_tmp.chromStart)[j])
            chromEnd.append(list(DF_tmp.chromEnd)[j])
        print 'finish',i
    data = {'chrom':chrom,'chromStart':chromStart,'chromEnd':chromEnd}
    DF = pd.DataFrame(data,columns=['chrom','chromStart','chromEnd'])
    return DF

def paixu2(DF):
    a = ['chr'+str(n) for n in range(1,23)]
    a.append('chrX')
    a.append('chrY')
    chrom = [];
    chromStart = [];
    chromEnd = [];
    length = [];
    signalValue = [];
    EID = [];
    for i in a:
        DF_tmp = DF.loc[DF["chrom"] == i].sort('chromStart')
        print 'begin',i
        for j in range(0,len(list(DF_tmp.chrom)),1):
            chrom.append(list(DF_tmp.chrom)[j])
            chromStart.append(list(DF_tmp.chromStart)[j])
            chromEnd.append(list(DF_tmp.chromEnd)[j])
            length.append(list(DF_tmp.length)[j])
            signalValue.append(list(DF_tmp.signalValue)[j])
            EID.append(list(DF_tmp.EID)[j])
        print 'finish',i
    data = {'chrom':chrom,'chromStart':chromStart,'chromEnd':chromEnd,'length':length,
            'signalValue':signalValue,'EID':EID}
    DF = pd.DataFrame(data,columns=['chrom','chromStart','chromEnd','length','signalValue','EID'])
    return DF

def paixu3(DF):
    a = ['chr'+str(n) for n in range(1,23)]
    a.append('chrX')
    a.append('chrY')
    chrom = [];
    TSS = [];
    gene_name = [];
    strand = [];
    for i in a:
        DF_tmp = DF.loc[DF["chrom"] == i].sort('TSS')
        print 'begin',i
        for j in range(0,len(list(DF_tmp.chrom)),1):
            chrom.append(list(DF_tmp.chrom)[j])
            TSS.append(list(DF_tmp.TSS)[j])
            gene_name.append(list(DF_tmp.gene_name)[j])
            strand.append(list(DF_tmp.strand)[j])
        print 'finish',i
    data = {'chrom':chrom,'TSS':TSS,'gene_name':gene_name,'strand':strand,}
    DF = pd.DataFrame(data,columns=['chrom','TSS','gene_name','strand'])
    return DF

def get_enhancer():
    path = os.path.join(get_data_dir(), "tmp", "enhancer.txt");
    handle = SeqIO.parse(path,'fasta')
    chrom = [];
    chromStart = [];
    chromEnd = [];
    for i in handle:
        chrom.append(re.split(':|-|\\|',i.id)[1])
        chromStart.append(int(re.split(':|-|\\|',i.id)[2]))
        chromEnd.append(int(re.split(':|-|\\|',i.id)[3]))
    data = {'chrom':chrom,'chromStart':chromStart,'chromEnd':chromEnd}
    enh = pd.DataFrame(data,columns=['chrom','chromStart','chromEnd'])
    enh = paixu(enh)
    enh['chromMiddle'] = (enh.chromStart + enh.chromEnd)/2
    enh.to_csv(os.path.join(get_data_dir(), "tmp", "enh.csv"), 
              sep='\t', index=False)
              
def get_TSS():
    '''
    The predicted TSS was download at http://fantom.gsc.riken.jp/5/datafiles/phase1.3/extra/TSS_classifier/.
    '''
    path = os.path.join(get_data_dir(), "tmp", "TSS_human.bed");
    TSS = pd.read_csv(path, sep='\t',header=None)
    TSS.convert_objects(convert_numeric=True)
    # print type(TSS[1])
    data = {'chrom':list(TSS[0]),'chromStart':list(TSS[1]),'chromEnd':list(TSS[2])}
    TSS = pd.DataFrame(data,columns=['chrom','chromStart','chromEnd'])
    # TSS = paixu(TSS)
    TSS['chromMiddle'] = (TSS.chromStart + TSS.chromEnd)/2
    print TSS
    '''
    TSS.to_csv(os.path.join(get_data_dir(), "tmp", "TSS.csv"), 
              sep='\t', index=False)
              '''
def get_TSS2(): 
    path = os.path.join(get_data_dir(), "tmp", "Ensembl_v65.Gencode_v10.ENSG.csv");
    TSS = pd.read_csv(path, sep=',',header=None)
    TSS = TSS.loc[TSS[5] == 'protein_coding']
    chrom = list('chr'+TSS[1])
    chromStart = list(TSS[2])
    chromEnd = list(TSS[3])
    gene_name = list(TSS[6])
    positive = list(TSS[4])
    Tss = []
    for i in range(0,len(TSS.index),1):
        if positive[i] == 1:
            Tss.append(chromStart[i])
        else:
            Tss.append(chromEnd[i])
    data = {'chrom':chrom,'TSS':Tss,'gene_name':gene_name,'strand':positive}
    TSS = pd.DataFrame(data,columns=['chrom','TSS','gene_name','strand'])
    TSS.to_csv(os.path.join(get_data_dir(), "tmp", "TSS.csv"), 
              sep='\t', index=False)


def map_mark_state(mark,state,cut_off,TSS_file='TSS.csv',gene_id="gene_name"):
    a = ['chr'+str(n) for n in range(1,23)]
    a.append('chrX')
    a.append('chrY')
    TSS = os.path.join(get_data_dir(), "tmp", TSS_file)
    TSS_DF = pd.read_csv(TSS, sep='\t')
    
    EID = []
    chrom = []
    Tss = []
    signalValue = []
    chromStart = []
    chromEnd = []
    gene_name = []
    '''
    path = os.path.join(get_data_dir(), "tmp", "{0}-{1}.csv".format('chr1', mark))
    DF = pd.read_csv(path, sep='\t')
    TSS_DF_tmp = TSS_DF.loc[TSS_DF['chrom'] == 'chr1']
    middle = list((DF['chromStart'].values+DF['chromEnd'].values)/2)
    print middle
    '''
    for c in a:
        path = os.path.join(get_data_dir(), "tmp", "{0}-{1}.csv".format(c, mark))
        DF = pd.read_csv(path, sep='\t')
        TSS_DF_tmp = TSS_DF.loc[TSS_DF['chrom'] == c]
        middle = list((DF['chromStart'].values+DF['chromEnd'].values)/2)
        DF['middle'] = middle
        for i,k in zip(list(TSS_DF_tmp['TSS'].values),list(TSS_DF_tmp[gene_id])):
            tmp = DF.loc[i-1000<DF['middle']]
            tmp = tmp.loc[DF['middle']<i+1000] 
            if len(tmp.index) >= cut_off :
                for j in range(0,len(tmp.index),1):
                    chrom.append(c)
                    Tss.append(i)
                    EID.append(list(tmp.EID)[j])
                    signalValue.append(list(tmp.signalValue)[j])
                    chromStart.append(list(tmp.chromStart)[j])
                    chromEnd.append(list(tmp.chromEnd)[j])
                    gene_name.append(k)
                
        print 'finish',c
               
    data = {'chrom':chrom,'chromStart':chromStart,'chromEnd':chromEnd,'TSS':Tss,'EID':EID,'signalValue':signalValue,'gene_name':gene_name}
    DF = pd.DataFrame(data,columns=['chrom','chromStart','chromEnd','TSS','EID','signalValue','gene_name'])
    DF.to_csv(os.path.join(get_data_dir(), "tmp", "{0} in {1}-{2}.csv".format(mark,state,cut_off)),
               sep='\t', index=False)
           
def he_bing_feng(mark,state,cut_off):
    path1 = os.path.join(get_data_dir(), "tmp", "{0} in {1}-{2}.csv".format(mark,state,cut_off))
    DF1 = pd.read_csv(path1, sep='\t')
    path2 = os.path.join(get_data_dir(), "tmp", 'TSS.csv')
    DF2 = pd.read_csv(path2,sep='\t')
    DF2 = paixu3(DF2)
    Full_EID_list = get_full_EID_list()
    width = [];
    signalValue = [];
    chrom = [];
    TSS = []
    EID = []
    gene_name = []
    n = 0
    
    for i,j,k in zip(list(DF2['chrom']),list(DF2['TSS']),list(DF2['gene_name'])):
        print n 
        n += 1
        tmp = DF1.loc[(DF1['chrom'] == i) & (DF1['TSS'] == j)]
        for id in Full_EID_list:
            tmp2 = tmp.loc[ tmp['EID'] == id] 
            if len(tmp2.index) >= 2:
                signalValue.append(sum(tmp2['signalValue'].values)/len(tmp2.index))
                width.append((sum(tmp2['chromEnd'].values)-sum(tmp2['chromStart'].values))/len(tmp2.index))
            elif len(tmp2.index) == 1:
                signalValue.append(list(tmp2['signalValue'])[0])
                width.append(list(tmp2['chromEnd'])[0]-list(tmp2['chromStart'])[0])
            else:
                continue
            chrom.append(i)
            TSS.append(j)
            EID.append(id)
            gene_name.append(k)
            
    data = {'chrom':chrom,'TSS':TSS,'width':width,'EID':EID,'signalValue':signalValue,'gene_name':gene_name}
    DF = pd.DataFrame(data,columns=['chrom','TSS','width','EID','signalValue','gene_name'])
    DF.to_csv(os.path.join(get_data_dir(), "tmp", "{0} in {1}-{2}organized.csv".format(mark,state,cut_off)),
               sep='\t', index=False)
    
def res_matrix(mark,state,cut_off=40):
    path = os.path.join(get_data_dir(), "tmp", "{0} in {1}-{2}.csv".format(mark, state,cut_off))
    DF = pd.read_csv(path, sep='\t')
    Full_EID_list = get_full_EID_list()
    res_matrix = []
    tmp = [0.]*len(Full_EID_list)
    for i in range(0,len(DF.index),1):
        try:
            if DF.chromMiddle[i-1] == DF.chromMiddle[i]:
                tmp[Full_EID_list.index(DF.EID[i])] = DF.signalValue[i]
            else:
                res_matrix.append(tmp)
                tmp = [0.]*len(Full_EID_list)
        except:
            pass
    
    f, ax = plt.subplots(figsize=(15, 15))
    cmap = sns.diverging_palette(210, 10, as_cmap=True)
    sns.corrplot(np.array(res_matrix), annot=False, sig_stars=False,   # .T??
             diag_names=False, cmap=cmap, ax=ax)
    f.tight_layout()
    plt.show()
     
    path2 = os.path.join(get_data_dir(), "tmp","{0} in {1}-{2}_diff.csv".format(mark,state,cut_off)) 
    a = open(path2,'w')
    for i in range(0,len(res_matrix[0]),1):
        for j in range(0,len(res_matrix),1):
            a.write(str(res_matrix[j][i])+"\t")
        a.write("\n")
    a.close() 
    
def organized():
    path = os.path.join(get_data_dir(), "tmp", "top width -0.01csv");
    width = pd.read_csv(path, sep='\t')
    width = paixu2(width)
    width.to_csv(os.path.join(get_data_dir(), "tmp", 'top width -0.01_organized.csv'),
               sep='\t', index=False)
    
    
def tmp_justhaveatry():  
    path = os.path.join(get_data_dir(), "tmp", "H3K4me3 in TSS-40organized.csv");
    DF = pd.read_csv(path, sep='\t')
    DF = DF.sort('width',ascending=False)
    # DF = DF.loc[DF['EID'] == EID]
    x = DF['width'].values
    y = DF['signalValue'].values 
    print DF
    m = [(i+1.0)/1 for i in range(0,100,1)]
    n = [];
    chang = [];
    for i in m:
        n.append(sum(y[len(y)*(i-1)/100:len(y)*(i)/100])/(len(y)/100))
        chang.append((sum(x[len(x)*(i-1)/100:len(x)*(i)/100])/(len(x)/100)))
    data = {'width':chang,'signalValue':n}
    print chang
    DF = pd.DataFrame(data,columns=['width','signalValue'])
    DF.to_csv(os.path.join(get_data_dir(), "tmp", "width vs signal.csv"),
               sep='\t', index=False)
    plot(chang,n,'.')
    plt.show()
        
    
    
if __name__ == "__main__":
    # get_enhancer()
    # get_TSS2()
    mark = 'H3K4me3'
    state = 'TSS'
    cut_off = 40 
    map_mark_state(mark,state,cut_off,TSS_file='ENSG_TSS.csv',gene_id="ENSG_ID");
    # res_matrix(mark,state)
    # organized();
    # tmp_justhaveatry()
    # he_bing_feng(mark,state,cut_off)
