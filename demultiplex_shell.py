#! /usr/env python3

import Demultiplex
import time


def launch_demultiplex(*args, directory='path', sample_key='path', mismatch=1, file_label='rbbr', barcode_1=None,
                 barcode_2=None, output_directory=None):

    start_time = time.time()
    demultiplex = Demultiplex.Demuliplex(*args, directory=directory, barcode_1=barcode_1, barcode_2=barcode_2,
                                         sample_key=sample_key, file_label=file_label)
    demultiplex.get_directory_lists()
    demultiplex.process_barcodes()
    demultiplex.process_file_label()
    demultiplex.get_sample_labels()
    demultiplex.output_objects(output_directory=output_directory)
    demultiplex.iterate_through_gseq()
    end_time = time.time()
    
    print('Total readss:' + str(demultiplex.reads))
    print('Reads passing filter:' + str(demultiplex.reads_pass_filter))
    print('Unmatched reads:' + str(demultiplex.unmatched_read))
    print('Total time:' + str(round((end_time - start_time) / 60.0, 2)) + ' minutes')

launch_demultiplex('s_1_1_*_qseq.txt', 's_1_2_*_qseq.txt', 's_1_3_*_qseq.txt', 's_1_4_*_qseq.txt',
                      directory='/Users/colinfarrell/Desktop/demultiplex_test_folder/',
                      barcode_1='/Users/colinfarrell/Desktop/demultiplex_test_folder/N700_nextera_bacrodes.txt',
                      barcode_2='/Users/colinfarrell/Desktop/demultiplex_test_folder/nextera_n50x_barcodes.txt',
                      sample_key='/Users/colinfarrell/Desktop/demultiplex_test_folder/atac_sample_key.txt',
                   output_directory='/Users/colinfarrell/Desktop/demultiplex_test_folder/output_test/')

# /u/scratch2/c/colinpat/ATACseq/demultiplexed
