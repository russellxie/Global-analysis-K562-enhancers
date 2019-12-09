#!/usr/bin/env python3
import os
import sys
import collections
import argparse
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

from _preprocessing import *
from _util import *

np.random.seed(0)

print("Start Analyzing.", file = sys.stderr)

def get_matrix_from_h5(filename, genome):
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

def save_matrix_to_h5(gbm, filename, genome):
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

def subsample_matrix(gbm, barcode_indices):
    return GeneBCMatrix(gbm.gene_ids,
                        gbm.gene_names,
                        gbm.barcodes[barcode_indices],
                        gbm.matrix[:, barcode_indices])

def get_expression(gbm, gene_name):
    gene_indices = np.where(gbm.gene_names == gene_name)[0]
    if len(gene_indices) == 0:
        raise Exception("%s was not found in list of gene names." % gene_name)
    return gbm.matrix[gene_indices[0], :].toarray().squeeze()

#Create an N-dimentional nested dictionary
def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

#Function to define the turning point of a cumulative curve by using the minimum derivative method
def turn_point(sgRNA_name, df):
    sgRNA_count  = df.T.filter(items=[sgRNA_name]).sum(axis=1).sort_values(ascending=False)
    sgRNA_cumsum = sgRNA_count.cumsum()

    #get the total cell number of this sgRNA
    cell_num = np.argwhere(sgRNA_count > 0).size

    #calculate the turning point by using the max derivative
    turning_point = sgRNA_cumsum.loc[((sgRNA_cumsum.diff()) / sgRNA_count.sum() > (1/cell_num))].shape

    return(sgRNA_count.iloc[turning_point])

#load sgRNA-barcode data
def load_data(input_file):
    input_fh  = open(input_file, 'r')
    #define parameters
    cell_bc_list   = []
    num_sgRNA_list = np.array([])
    sgRNAs         = []
    umis           = []
    #initiate a 2D dictionary
    data_dict = nested_dict(2, list)
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

#filter the sgRNA UMI count based on the cutoff values
def filter_umi (df):
    #calculate the cutoff value for each sgRNA in the dataset
    sgRNA_cutoff = [[turn_point(i, df)] for i in list(df.index)]
    for i in range(0, len(sgRNA_cutoff)):
        df.iloc[i, :].loc[df.iloc[i, :] <= sgRNA_cutoff[i]] = 0
    return df

#fine the indecies of cells which contain the input sgRNA
def find_sgrna_cells(sgRNA_list, df, gbm):
    cell_bc = df.loc[:, (df.loc[sgRNA_list].sum() > 0)].T.index
    cell_index = []
    for i in cell_bc:
        current_idx = np.argwhere(gbm.barcodes == i.encode('utf-8'))
        if  current_idx.size > 0:
            cell_index.append(current_idx.item())
    return [x for x in set(cell_index)]

def find_non_zero_cells(df, gbm):
    cell_bc = df.T.loc[df.sum() != 0].index
    cell_index = []
    for i in cell_bc:
        current_idx = np.argwhere(gbm.barcodes == i.encode('utf-8'))
        if  current_idx.size > 0:
            cell_index.append(current_idx.item())
    return [x for x in set(cell_index)]

def hypergeo_test(non_zero_array, sgrna_idx, i):
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
def perform_DE(sgrna_idx, input_array, idx, num_processes, pval_list_down, pval_list_up, fc_list):
    nonzero_pval_list_up = []
    nonzero_pval_list_down = []
    nonzero_fc_list = []
    with Pool(processes=num_processes) as p:
        for pval_down, pval_up, fc in p.starmap(hypergeo_test, zip(
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

def FDR(x):
    """
    Assumes a list or numpy array x which contains p-values for multiple tests
    Copied from p.adjust function from R  
    """
    o = [i[0] for i in sorted(enumerate(x), key=lambda v:v[1],reverse=True)]
    ro = [i[0] for i in sorted(enumerate(o), key=lambda v:v[1])]
    q = sum([1.0/i for i in range(1,len(x)+1)])
    l = [q*len(x)/i*x[j] for i,j in zip(reversed(range(1,len(x)+1)),o)]
    l = [l[k] if l[k] < 1.0 else 1.0 for k in ro]
    return np.asarray(l)

def create_combo(sgrna_list):
    combo_list = []
    for i in range(1, 11):
        my_array = list(itertools.combinations(enumerate(sgrna_list), i))
        for j in range(len(my_array)):
            idx = list(zip(*my_array[j]))[0]
            combo_list.append(idx)
    #return a list of tuple, which idicate the indexing of each combination
    return combo_list;

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', '--hdf5', dest='input_hdf5', required=True,
        type=str,
        help='specify the hdf5 file (output from 10X pipeline).'
    )
    parser.add_argument(
        '-s', '--sgrna', dest='input_sgrna', required=True,
        type=str,
        help='specify the sgrna summary file.'
    )
    parser.add_argument(
        '-g', '--genome', dest='input_genome', required=True,
        type=str,
        help='specify the genome of the hdf5 file.'
    )
    parser.add_argument(
        '-l', '--sgrna_list', dest='sgrna_list', required=True,
        type=str,
        help='specify the sgRNAs need to be tested.'
    )
    parser.add_argument(
        '-o', '--output_dir', dest='output_dir', required=False,
        default = '.',
        help='specify an output directory, default is current dir.'
    )
    parser.add_argument(
        '-t', '--threads', dest = 'threads', required=False,
        type=int,
        default=1,
        help='set number of barcode comparison threads. \
                        The default is 1'
    )
    parser.add_argument(
        '-n', '--norm', dest = 'norm_method', required=False,
        type=str,
        default='cpm',
        help='choose normalization methods: CPM only or \
                        normalize to the Meta-cell.'
    )

    args = parser.parse_args()

    num_processing = args.threads
    norm = args.norm_method

    #check the normalization method
    if (norm != 'cpm') and (norm != 'metacell'):
        print("Incorrect normalization method. Has to be either 'cpm' or 'metacell'.", file = sys.stderr)
        sys.exit(0)

    #read the 10X hdf5 file
    print("Reading HDF5 File...", file = sys.stderr)
    filtered_matrix_h5 = args.input_hdf5
    genome = args.input_genome
    gene_bc_matrix = get_matrix_from_h5(filtered_matrix_h5, genome)

    #load the sgRNA file
    print("Loading sgRNA Annotations...", file = sys.stderr)
    sgrna_file = args.input_sgrna
    sgrna_df_adj = pd.read_pickle(sgrna_file)

    #check the positive control sgRNA in FADS1 region, which has 287 cells
    #perform hypergeometric test for every single gene in the dataframe
    sgrna_dict  = {}
    sgrnas_file = args.sgrna_list
    with open(sgrnas_file) as f:
        for line in f:
            region_id, sgrna_string = line.strip().split("\t")
            sgrnas = sgrna_string.split(";")
            sgrna_dict.update({region_id : sgrnas})

    output_dir = args.output_dir

    #get rid of the cells without any detectable sgRNA
    print("Removing Cells without sgRNA...", file = sys.stderr)
    non_zero_cell_idx = find_non_zero_cells (sgrna_df_adj, gene_bc_matrix)
    nonzero_gbm = subsample_matrix(gene_bc_matrix, non_zero_cell_idx)

    '''
    Read parameters for figure plotting
    '''
    print("Loading Plotting Parameters...", file = sys.stderr)
    cmap = matplotlib.cm.get_cmap('Set2')

    #create cumsum list of the chromsomal position
    length_list = [
        0, 248956422,491149951,689445510,879660065,1061198324,
        1232004303,1391350276,1536488912,1674883629,1808681051,
        1943767673,2077042982,2191407310,2298451028,2400442217,
        2490780562,2574038003,2654411288,2713028904,2777473071,
        2824183054,2875001522
    ]

    chr_order = [
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
        'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
        'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
        'chr21', 'chr22', 'chrX'
    ]

    #read the gene annotation file for plotting, which re-ordered the genes based on their
    #relative positions on the genome
    annot_df = pd.read_csv(
        './generate_annotations/plot_annotation.txt',
        header = None,
        sep='\t',
        names = ["idx", "gene_names", "chromosome", "pos", "strand", "color_idx", "chr_idx"]
    )

    #Normalization.
    [g, c] = nonzero_gbm.matrix.shape
    cpm_matrix = np.zeros((g, c))
    uniq_id = set()

    for x in nonzero_gbm.barcodes:
        uniq_id.add(x.decode()[-1])

    if (norm == 'cpm'):
        cpm_matrix = np.array(nonzero_gbm.matrix / nonzero_gbm.matrix.sum(axis = 0) * 1e6)
    elif (norm == 'metacell'):
        for lib in sorted(uniq_id):
            print("Normalizing Batch " + lib, file = sys.stderr)
            index = [i for i,e in enumerate(nonzero_gbm.barcodes) if e.decode()[-1] == lib]
            one_cell_cpm = np.array((nonzero_gbm.matrix[:,index].sum(axis = 1) + 1)
                                    / np.sum(nonzero_gbm.matrix[:,index]) * 1e6).flatten()
            for i in index:
                cell_cpm = nonzero_gbm.matrix[:,i].toarray().flatten()\
                           / np.sum(nonzero_gbm.matrix[:,i]) * 1e6
                cpm_matrix[:,i] = cell_cpm / one_cell_cpm

    #prepare the input array
    nonzero_idx = np.where(np.sum(nonzero_gbm.matrix > 0, axis = 1) > 1)[0]
    idx = list(set(nonzero_idx) & set(annot_df.idx))

    #create input ndarray
    input_array = cpm_matrix[idx]

    for k in sgrna_dict:

        print("Region in processing: " + k[0:], file = sys.stderr)

        #idx index of cells containing the given sgRNA
        sgrna_idx = find_sgrna_cells(sgrna_dict[k], sgrna_df_adj, nonzero_gbm)

        #force the up-tail p-vals of all zero expressed genes to be zero. (actual p-val is 1)
        pval_list_down = np.zeros(len(nonzero_gbm.gene_names))
        pval_list_up = np.zeros(len(nonzero_gbm.gene_names))

        fc_list = np.ones(len(nonzero_gbm.gene_names))

        #perform the differential gene analysis by using Virtual FACS
        num_sgrna_cell, pval_list_up, pval_list_down, fc_list = perform_DE(
            sgrna_idx,
            input_array,
            idx,
            num_processing,
            pval_list_down,
            pval_list_up,
            fc_list
        )

        #save all the output
        io.savemat(
            '%s/%s-up_log-pval'%(output_dir, k[0:]),
            {'matrix':pval_list_up}
        )

        io.savemat(
            '%s/%s-down_log-pval'%(output_dir, k[0:]),
            {'matrix':pval_list_down}
        )


if __name__ == '__main__':
    main()
