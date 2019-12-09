#!/usr/bin/env python3
'''
This script will take the paired-end sequencing results for the sgRNA enrichment library
and extract the cell barcodes and the sgRNA sequences by comparing with the
known sgRNA/barcode list. The input file are the two fastq files of each read.
'''
import argparse
import sys
import regex
import itertools
import numpy as np
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from multiprocessing import Pool, Value

__author__  = "Jialei Duan and Russell Xie"
__version__ = "0.0.2"

def get_read_sequence(fastq_iter1, fastq_iter2):
    for (_x1,y1,_z1), (_x2,y2,_z2) in zip(fastq_iter1, fastq_iter2):
        yield y1,y2

def compile_regex_known_list(known_list, num_mismatches):
    known_list_regex = [regex.compile('(?e)' + 
                                      '(' + i.rstrip() + ')' +
                                      '{e<=' + str(num_mismatches) + '}') for i in known_list]
    return known_list_regex

def compare_with_known_list(obs_sequences, known_barcodes_regex, known_sgRNAs_regex, known_barcodes, known_sgRNAs):

    barcode_match = list_compare_helper_bc(obs_sequences[0], known_barcodes, known_barcodes_regex)
    sgRNA_match   = list_compare_helper_sp(obs_sequences[1], known_sgRNAs, known_sgRNAs_regex)
    
    return([barcode_match, sgRNA_match])

def list_compare_helper_bc(sequence, known_list, known_list_regex):
    matched_info = []
    for i in known_list:
        if i == sequence[0:16]:
            matched_info = [sequence, [i], [(0,16)], [(0,0,0)]]
            break
        
    if len(matched_info) == 0:
        matched = [(i.search(sequence), i)
                   for i in known_list_regex if i.search(sequence)]

        matched_info = [sequence,
                        [i[1].pattern[i[1].pattern.rfind("(") + 1 :
                         i[1].pattern.rfind(")")] for i in matched],
                        [i[0].span() for i in matched],
                        [i[0].fuzzy_counts for i in matched]]

    return (matched_info)

def list_compare_helper_sp(sequence, known_list, known_list_regex):
    matched_info = []
    for i in known_list:

        if i in sequence:
            matched_info = [sequence, [i], [(19,38)], [(0,0,0)]]
            break
    
    if len(matched_info) == 0:
        matched = [(i.search(sequence), i)
                   for i in known_list_regex if i.search(sequence)]

        matched_info = [sequence,
                        [i[1].pattern[i[1].pattern.rfind("(") + 1 :
                         i[1].pattern.rfind(")")] for i in matched],
                        [i[0].span() for i in matched],
                        [i[0].fuzzy_counts for i in matched]]

    return (matched_info)


def format_output(i):
    i = list(i)
        
    if i[1]:
        i[1:4] = [','.join(i[1]),
                  ','.join(':'.join([str(iii) for iii in ii]) for ii in i[2]),
                  ','.join(':'.join([str(iii) for iii in ii]) for ii in i[3])]
    else:
        i[1:4] = ['NA'] * 3
    return i        

def get_fastq_file_handle(fq_file):
    if fq_file.endswith('gz'):
        import gzip
        handle = gzip.open(fq_file, mode='rt')
    elif fq_file.endswith('bz2'):
        import bz2
        handle = bz2.open(fq_file, mode='rt')
    else:
        handle = open(fq_file, 'r')
    return handle

#make multithreading
def make_multithread(inner_func, numthreads):
    def func_mt(*args):
        length = len(args[0])
        result = np.empty(length, dtype=np.float64)
        args = (result,) + args        
        chunklen = (length + 1) // numthreads
        chunks = [[arg[i * chunklen:(i + 1) * chunklen] for arg in args]
                  for i in range(numthreads)]
        ''' 
        You should make sure inner_func is compiled at this point, because
        the compilation must happen on the main thread. This is the case
        in this example because we use jit().
        '''
        threads = [threading.Thread(target=inner_func, args=chunk)
                   for chunk in chunks[:-1]]
        for thread in threads:
            thread.start()

        # the main thread handles the last chunk
        inner_func(*chunks[-1])

        for thread in threads:
            thread.join()
        return result
    return func_mt

def main():
    
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--input_fastq', dest='input_fastq',
                        required=True, nargs='*',
                        type=str,
                        help='specify a (paired) fastq file(s) as input')

    parser.add_argument('-b', '--input_barcodes', dest='input_barcodes',
                        required=True,
                        type=str,
                        help='specify a reference barcodes file as input')
    
    parser.add_argument('-r', '--input_sgrnas', dest='input_sgrnas',
                        required=True,
                        type=str,
                        help='specify a file containing known sgRNAs as input')

    parser.add_argument('-m', '--num_mismatches', required=False,
                        type=int,
                        default=1,
                        help='specify a fuzzy matching search threshold. \
                        The default is 1')

    parser.add_argument('-t', '--threads', required=False,
                        type=int,
                        default=1,
                        help='set number of barcode comparison threads. \
                        The default is 1')

    parser.add_argument('-o', '--output', dest='output', required=False,
                        help='specify an output file')

    args = parser.parse_args()


    if args.output and args.output != '-':
        out_fh = open(args.output, 'w')

    fq_files       = args.input_fastq
    barcodes_file  = args.input_barcodes
    sgrnas_file    = args.input_sgrnas
    num_processes  = args.threads
    num_mismatches = args.num_mismatches
        
    #check if there are two fastq files
    if len(fq_files) != 2:
        exit("Two fastq files required.\n")
        
    #read the fastq files
    fastq_iter1 = FastqGeneralIterator(get_fastq_file_handle(fq_files[0]))
    fastq_iter2 = FastqGeneralIterator(get_fastq_file_handle(fq_files[1]))

    #read known barcode files
    if barcodes_file.endswith('gz'):
        import gzip
        bc_handle = gzip.open(barcodes_file, mode='rt')
        known_barcodes = [i.rstrip().replace('-1', '') for i in bc_handle]
    else:
        known_barcodes = [i.rstrip().replace('-1','') for i in open(barcodes_file, 'r')]
    known_barcodes_regex = compile_regex_known_list(known_barcodes, num_mismatches)

    #read known sgRNA files
    known_sgRNAs = [i.rstrip() for i in open(sgrnas_file, 'r')]
    known_sgRNAs_regex = compile_regex_known_list(known_sgRNAs, num_mismatches)

    counter = Value('i', 0)
    if num_processes > 1:
        
        chunk_size = 10000
        items = list(itertools.islice(get_read_sequence(fastq_iter1, fastq_iter2), chunk_size))
            
        with Pool(processes = num_processes) as p:
            
            while items:
                counter.value = counter.value + 1000
                if ((counter.value % 1000 == 0) and (counter.value % 10000 != 0)):
                    print(".", end='', flush=True, file=sys.stdout)
                else:
                    print(counter.value, flush=True, file=sys.stdout)

                
                results = p.starmap(compare_with_known_list, zip(items,
                                                                 itertools.repeat(known_barcodes_regex),
                                                                 itertools.repeat(known_sgRNAs_regex),
                                                                 itertools.repeat(known_barcodes),
                                                                 itertools.repeat(known_sgRNAs)))
                for i in results:
                    if ((i[0][1]) and (i[1][1])):
                        print ('\t'.join(format_output(i[0])) + 
                               '\t' + 
                               '\t'.join(format_output(i[1])), file=out_fh)
            
                items = list(itertools.islice(get_read_sequence(fastq_iter1, fastq_iter2), chunk_size))

    else:
        for i in get_read_sequence(fastq_iter1, fastq_iter2):
            result = compare_with_known_list(i, known_barcodes_regex, known_sgRNAs_regex,
                                             known_barcodes, known_sgRNAs)

            if ((result[0][1]) and (result[1][1])):
                        print ('\t'.join(format_output(result[0])) + 
                               '\t' + 
                               '\t'.join(format_output(result[1])), file=out_fh)

    if args.output and args.output != '-':
        out_fh.close()

if __name__ == '__main__':
    main()
