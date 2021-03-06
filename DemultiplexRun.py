#!/usr/bin/env python3

import DemultiplexClass
import time
import argparse


def launch_demultiplex(*args, directory='path', sample_key='path', mismatch=1, file_label='', barcode_1=None,
                       barcode_2=None, output_directory=None, gnu_zipped=False):
    """Simple function to initialize DemultiplexClass"""
    start_time = time.time()
    demultiplex = DemultiplexClass.Demuliplex(*args, directory=directory, barcode_1=barcode_1, barcode_2=barcode_2,
                                              sample_key=sample_key, file_label=file_label, mismatch=mismatch,
                                              gnu_zipped=gnu_zipped)
    demultiplex.get_directory_lists()
    demultiplex.process_barcodes()
    demultiplex.process_file_label()
    demultiplex.get_sample_labels()
    demultiplex.output_objects(output_directory=output_directory)
    demultiplex.iterate_through_qseq()
    end_time = time.time()
    print('Total reads:' + str(demultiplex.reads))
    print('Reads passing filter:' + str(demultiplex.reads_pass_filter))
    print('Indexed reads:' + str(demultiplex.indexed_reads))
    print('Unmatched reads:' + str(demultiplex.unmatched_read))
    print('Total time:' + str(round((end_time - start_time) / 60.0, 2)) + ' minutes')


parser = argparse.ArgumentParser(description='Demultiplexing script. Script demultiplexes Illumnina qseq lane files '
                                             'outputing sample fastq files. Works with .gz and uncompressed qseq files.'
                                             ' Options for single and dual indexes'
                                             '\n\n Usage; demultiplex -D directory -S sample_key'
                                             ' -B1 barcode_1 -B2 barcode_2 -L file_labels -M mismatch_number -O '
                                             'output_directory -I input_file_1 input_file_2 ...')

parser.add_argument('-D', type=str, help='/path/ to qseq directory')
parser.add_argument('-S', type=str, help='/path/sample_file.txt file should be formatted as \''
                                         'barcode tab sample_name\' for single index and '
                                         '\'barcode tab barcode tab sample_name\' '
                                         'for dual indexes ')
parser.add_argument('-B1', type=str, help='/path/barcode_1_file, barcode \t index key')
parser.add_argument('-B2', type=str, default=None, help='/path/barcode_2_file, barcode \t index key')
parser.add_argument('-L', type=str, help='string of r and b character to designate input files as '
                                         'barcode or read files, should be the same order as input'
                                         'file')
parser.add_argument('-M', type=int, default=2, help='number of barcode mismatches to consider, default = 2')
parser.add_argument('-O', type=str, help='path to output directory')
parser.add_argument('-Z', action="store_true", default=False, help='if qseq files gzipped, slows processing')
parser.add_argument('-I', type=str, nargs='*', help='qseq file prefix and suffix separated'
                                                    'by ^, ie. -I s_1_^.qseq.txt '
                                                    's_2_^.qseq.txt ')
arguments = parser.parse_args()

print('Started Job')

launch_demultiplex(*arguments.I, directory=arguments.D, barcode_1=arguments.B1, barcode_2=arguments.B2,
                   sample_key=arguments.S, output_directory=arguments.O, mismatch=arguments.M, gnu_zipped=arguments.Z,
                   file_label=arguments.L)
