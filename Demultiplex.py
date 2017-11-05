#!/usr/bin/env python3
import os
import sys
import time
from BarcodeFileParser import BarcodeFileParser
from IncrementingFileExtractor import IncrementingFileExtractor
from KeyFileParser import KeyFileParser
from QseqFileParser import QseqFileParser

class Demuliplex:

    def __init__(self, 
        files = [],
        sample_key_path = None, 
        primary_barcodes_path = None,
        secondary_barcodes_path = None,
        mismatch = 1):

        self.start_time = time.time()
        self.end_time = time.time()

        self.file_extractor = IncrementingFileExtractor(files)
        self.primary_barcode_parser = BarcodeFileParser(path = primary_barcodes_path,
                mismatch = mismatch)
        self.secondary_barcode_parser = BarcodeFileParser(path = secondary_barcodes_path,
                mismatch = mismatch)
        self.key_file_parser = KeyFileParser(path = sample_key_path)

    def run(self, 
        output_directory = os.getcwd(), 
        gnu_zipped = False):

        self.start_time = time.time()

        files = self.file_extractor.get_files()
        primary_barcodes = self.primary_barcode_parser.get_barcodes()
        secondary_barcodes = self.secondary_barcode_parser.get_barcodes()

        if secondary_barcodes:
            self.key_file_parser.set_parse_type('dual')

        sample_list = self.key_file_parser.get_sample_list()

        self.qseq_file_parser = QseqFileParser(barcode_list = [primary_barcodes, secondary_barcodes],
            output_directory = output_directory,
            sample_list = sample_list,
            read_count = self.file_extractor.get_read_count())
        
        self.qseq_file_parser.run(files = files, gnu_zipped = gnu_zipped)    

        self.end_time = time.time()

    def print_output(self):
        print('Total reads:' + str(self.qseq_file_parser.reads))
        print('Reads passing:' + str(self.qseq_file_parser.reads_pass))
        print('Indexed reads:' + str(self.qseq_file_parser.indexed_reads))
        print('Unmatched reads:' + str(self.qseq_file_parser.unmatched_reads))
        print('Total time:' + str(round((self.end_time - self.start_time) / 60.0, 2)) + ' minutes')