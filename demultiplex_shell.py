#! /usr/env python3

import Demultiplex
import time

def launch_demultiplex(*args, 
    directory = 'path', 
    sample_key = 'path', 
    mismatch = 1, 
    file_label = 'rbbr', 
    barcode_1 = None,
    barcode_2 = None, 
    output_directory = None):

    start_time = time.time()

    demultiplex = Demultiplex.Demuliplex(*args, 
        directory = directory, 
        output_directory = output_directory,
        barcode_1 = barcode_1, 
        barcode_2 = barcode_2,
        sample_key = sample_key, 
        file_label = file_label)

    demultiplex.run()
    end_time = time.time()
    demultiplex.print()

launch_demultiplex('s_1_1_*_qseq.txt', 
    's_1_2_*_qseq.txt', 
    's_1_3_*_qseq.txt', 
    's_1_4_*_qseq.txt',
    directory = '/Users/colinfarrell/Desktop/demultiplex_test_folder/',
    barcode_1 = '/Users/colinfarrell/Desktop/demultiplex_test_folder/N700_nextera_bacrodes.txt',
    barcode_2 = '/Users/colinfarrell/Desktop/demultiplex_test_folder/nextera_n50x_barcodes.txt',
    sample_key = '/Users/colinfarrell/Desktop/demultiplex_test_folder/atac_sample_key.txt',
    output_directory = '/Users/colinfarrell/Desktop/demultiplex_test_folder/output_test/')

# /u/scratch2/c/colinpat/ATACseq/demultiplexed
