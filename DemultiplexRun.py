#!/usr/bin/env python3

import Demultiplex
import time
import argparse

parser = argparse.ArgumentParser(description = 'Demultiplexing script. Script demultiplexes Illumnina qseq lane files '
                                             'outputing sample fastq files. Works with .gz and uncompressed qseq files.'
                                             ' Options for single and dual indexes'
                                             '\n\n Usage; demultiplex -D directory -S sample_key'
                                             ' -B1 barcode_1 -B2 barcode_2 -L file_labels -M mismatch_number -O '
                                             'output_directory -I input_file_1 input_file_2 ...')

parser.add_argument('-D', type = str, help = '/path/ to qseq directory')
parser.add_argument('-S', type = str, help = '/path/sample_file.txt file should be formatted as \''
                                         'barcode tab sample_name\' for single index and '
                                         '\'barcode tab barcode tab sample_name\' '
                                         'for dual indexes ')
parser.add_argument('-B1', type = str, help = '/path/barcode_1_file, barcode \t index key')
parser.add_argument('-B2', type = str, default = None, help = '/path/barcode_2_file, barcode \t index key')
parser.add_argument('-L', type = str, help = 'string of r and b character to designate input files as '
                                         'barcode or read files, should be the same order as input'
                                         'file')
parser.add_argument('-M', type = int, default = 2, help = 'number of barcode mismatches to consider, default = 2')
parser.add_argument('-O', type = str, help = 'path to output directory')
parser.add_argument('-Z', action = "store_true", default = False, help = 'if qseq files gzipped, slows processing')
parser.add_argument('-I', type = str, nargs = '*', help = 'qseq file prefix and suffix separated'
                                                    'by ^, ie. -I s_1_^.qseq.txt '
                                                    's_2_^.qseq.txt ')

arguments = parser.parse_args()

print('Started Job')

files = []

for index, file in enumerate(arguments.I):
    files.append({
        "path": arguments.D + '/' + file,
        "action": arguments.L[index]
    })

demultiplex = Demultiplex.Demuliplex(files = files,
        sample_key_path = arguments.S,
        primary_barcodes_path = arguments.B1,
        secondary_barcodes_path = arguments.B2,
        mismatch = arguments.M)

demultiplex.run(output_directory = arguments.O, gnu_zipped = arguments.Z)

demultiplex.print_output()