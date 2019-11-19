#!/usr/bin/env python3

import collections
import tables
import itertools

import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse

from multiprocessing import Pool
from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix
from multiprocessing import Pool

def _load_data(input_file):
    from . import _nested_dict

    '''
    Inputs
    ---------
    Load the sgRNA-barcode file for the correction

    The file should have the following columns:
    - 10X Cell-barcode
    - Total Reads
    - Total sgRNA Count
    - sgRNA sequqnce (seperated by colon)
    - Read count for each sgRNA (separated by colon)
    - UMI for each sgRNA (separated by colon)
                                   
    Example line:
    CGTAGGCGTTGGTTTG-2	2	2	CTGTTTTAGGACTTTAGAC;TTCCGCGTTACATAACTTA	1;1	1;1

    Returns
    --------
     dataframe which contains the matrix of all cell-barcodes vs all sgRNAs.

    '''
    input_fh  = open(input_file, 'r')
    #define parameters
    cell_bc_list   = []
    num_sgRNA_list = np.array([])
    sgRNAs         = []
    umis           = []
    #initiate a 2D dictionary
    data_dict = _nested_dict(2, list)
    for line in input_fh:
        cell_bc    = line.strip().split('\t')[0]
        num_sgRNA  = line.strip().split('\t')[2]
        sgRNA_list = line.strip().split('\t')[3].split(';')
        umi_list   = line.strip().split('\t')[5].split(';')
        for i in zip(sgRNA_list, umi_list):
            data_dict[cell_bc][i[0]] = i[1]
    #read the 2D dictionary into the pandas DataFrame
    df = pd.DataFrame(data_dict).fillna(0).astype(int)    
    return df

def _nested_dict(n, type):
    '''
    Create an N-dimentional nested dictionary
    '''

    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))


#filter the sgRNA UMI count based on the cutoff values
def _filter_umi (df):
    #calculate the cutoff value for each sgRNA in the dataset
    sgRNA_cutoff = [[turn_point(i, df)] for i in list(df.index)]
    for i in range(0, len(sgRNA_cutoff)):
        df.iloc[i, :].loc[df.iloc[i, :] <= sgRNA_cutoff[i]] = 0
    return df


#Function to define the turning point of a cumulative curve by using the minimum derivative method
def _turn_point(sgRNA_name, df):
    sgRNA_count  = df.T.filter(items=[sgRNA_name]).sum(axis=1).sort_values(ascending=False)
    sgRNA_cumsum = sgRNA_count.cumsum()

    #get the total cell number of this sgRNA
    cell_num = np.argwhere(sgRNA_count > 0).size

    #calculate the turning point by using the max derivative
    turning_point = sgRNA_cumsum.loc[((sgRNA_cumsum.diff()) / sgRNA_count.sum() > (1/cell_num))].shape
    
    return(sgRNA_count.iloc[turning_point])


