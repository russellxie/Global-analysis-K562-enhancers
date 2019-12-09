#!/usr/bin/env python3
import collections
import tables
import itertools

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import scipy.sparse as sp_sparse

from multiprocessing import Pool
from collections import defaultdict
from scipy import sparse, io
from scipy.sparse import csr_matrix
from multiprocessing import Pool

def _get_matrix_from_h5(filename, genome):
    '''
    From the 10X Cell ranger pipeline. 
    This code only works with the 10X V2 output
    '''
    GeneBCMatrix = collections.namedtuple('GeneBCMatrix',
                                      ['gene_ids', 'gene_names', 'barcodes', 'matrix'])

    with tables.open_file(filename, 'r') as f:
        try:
            dsets = {}
            for node in f.walk_nodes('/' + genome, 'Array'):
                dsets[node.name] = node.read()
            matrix = sp_sparse.csc_matrix((dsets['data'],
                                           dsets['indices'],
                                           dsets['indptr']),
                                           shape=dsets['shape'])
            return GeneBCMatrix(dsets['genes'],
                                dsets['gene_names'],
                                dsets['barcodes'],
                                matrix)
        except tables.NoSuchNodeError:
            raise Exception("Genome %s does not exist in this file." % genome)
        except KeyError:
            raise Exception("File is missing one or more required datasets.")

def _save_matrix_to_h5(gbm, filename, genome):
    flt = tables.Filters(complevel=1)
    with tables.open_file(filename, 'w', filters=flt) as f:
        try:
            group = f.create_group(f.root, genome)
            f.create_carray(group, 'genes', obj=gbm.gene_ids)
            f.create_carray(group, 'gene_names', obj=gbm.gene_names)
            f.create_carray(group, 'barcodes', obj=gbm.barcodes)
            f.create_carray(group, 'data', obj=gbm.matrix.data)
            f.create_carray(group, 'indices', obj=gbm.matrix.indices)
            f.create_carray(group, 'indptr', obj=gbm.matrix.indptr)
            f.create_carray(group, 'shape', obj=gbm.matrix.shape)
        except:
            raise Exception("Failed to write H5 file.")

def _get_expression(gbm, gene_name):
    gene_indices = np.where(gbm.gene_names == gene_name)[0]
    if len(gene_indices) == 0:
        raise Exception("%s was not found in list of gene names." % gene_name)
    return gbm.matrix[gene_indices[0], :].toarray().squeeze()


#fine the indecies of cells which contain the input sgRNA
def _find_sgrna_cells(sgRNA_list, df, gbm):
    cell_bc = df.loc[:, (df.loc[sgRNA_list].sum() > 0)].T.index
    cell_index = []
    for i in cell_bc:
        current_idx = np.argwhere(gbm.barcodes == i.encode('utf-8'))
        if  current_idx.size > 0:
            cell_index.append(current_idx.item())
    return [x for x in set(cell_index)]

def _hypergeom_test(non_zero_array, sgrna_idx, i):
    #find indecies of cells in which expression of given gene is
    #equal or less than the median of this gene in the whole population
    median_cell_idx  = np.argwhere(non_zero_array <= np.median(non_zero_array))

    #find the same cells subset in the cells with a given sgRNA
    overlap_cell_idx = np.intersect1d(median_cell_idx, sgrna_idx)

    #calculate the median fold change
    other_idx = np.setxor1d(sgrna_idx, range(len(non_zero_array)))

    fc = (np.mean(non_zero_array[sgrna_idx]) + 0.01) / (np.mean(non_zero_array[other_idx]) + 0.01)

    #perform hypergeometric test, get the upper tail
    k = len(overlap_cell_idx)
    M = len(non_zero_array)
    n = len(median_cell_idx)
    N = len(sgrna_idx)
    try:
        pval_up = stats.hypergeom.logcdf(k, M, n, N).item()
    except:
        pval_up = float('nan')

    try:
        pval_down = stats.hypergeom.logsf(k, M, n, N).item()
    except:
        pval_down = float('nan')

    return pval_down, pval_up, fc

#num_sgrna_cell, pval_list = perform_DE(sgrna_idx, input_array, idx, num_processing)
def _perform_DE(sgrna_idx, input_array, idx, num_processes, pval_list_down, pval_list_up, fc_list):
    '''
    Helper function to calculate the hypergeometric p value. 
    '''
    nonzero_pval_list_up = []
    nonzero_pval_list_down = []
    nonzero_fc_list = []
    with Pool(processes=num_processes) as p:
        for pval_down, pval_up, fc in p.starmap(_hypergeom_test, zip(
                input_array,
                itertools.repeat(sgrna_idx),
                idx)
        ):
            nonzero_pval_list_down.append(pval_down)
            nonzero_pval_list_up.append(pval_up)
            nonzero_fc_list.append(fc)
    for i in idx:
        pval_list_up[i] = nonzero_pval_list_up.pop(0)
        pval_list_down[i] = nonzero_pval_list_down.pop(0)
        fc_list[i] = nonzero_fc_list.pop(0)
    return len(sgrna_idx), pval_list_up, pval_list_down, fc_list

def _FDR(x):
    """
    -Deprecated, now we calculate the backgroud p-val for every gene explicitly.

    Benjamini-Hochberg FDR calculation.
    Assumes a list or numpy array x which contains p-values for multiple tests
    Copied from p.adjust function from R
    """
    o = [i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True)]
    ro = [i[0] for i in sorted(enumerate(o), key=lambda v:v[1])]
    q = sum([1.0/i for i in range(1,len(x)+1)])
    l = [q*len(x)/i*x[j] for i,j in zip(reversed(range(1,len(x)+1)),o)]
    l = [l[k] if l[k] < 1.0 else 1.0 for k in ro]
    return np.asarray(l)

def _create_combo(sgrna_list):
    combo_list = []
    for i in range(1, 11):
        my_array = list(itertools.combinations(enumerate(sgrna_list), i))
        for j in range(len(my_array)):
            idx = list(zip(*my_array[j]))[0]
            combo_list.append(idx)
    #return a list of tuple, which idicate the indexing of each combination
    return combo_list

