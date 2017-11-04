#!/usr/bin/env python3
import os
import sys
from BarcodeFileParser import BarcodeFileParser
from FileExtractor import FileExtractor
from KeyFileParser import KeyFileParser
from QseqFileParser import QseqFileParser


class Demuliplex:

    def __init__(self, 
        files = [],
        sample_key_path = None, 
        mismatch = 1, 
        primary_barcodes_path = None,
        secondary_barcodes_path = None):

        self.file_extractor = FileExtractor(files)
        self.primary_barcode_parser = BarcodeFileParser(path = primary_barcodes_path,
                mismatch = mismatch)
        self.secondary_barcode_parser = BarcodeFileParser(path = secondary_barcodes_path,
                mismatch = mismatch)
        self.key_file_parser = KeyFileParser(path = sample_key_path)

    def run(self, 
        output_directory = os.getcwd(), 
        gnu_zipped = False):

        files = self.file_extractor.get_files()
        primary_barcodes = self.primary_barcode_parser.get_barcodes()
        secondary_barcodes = self.secondary_barcode_parser.get_barcodes()

        parse_type = 'single'
        if secondary_barcodes:
            parse_type = 'dual'
            
        self.key_file_parser.set_labels(parse_type)

        sample_list = self.key_file_parser.get_sample_list()
        sample_key = self.key_file_parser.get_sample_key()

        qseq_file_parser = QseqFileParser(files = files,
            barcode_list = [primary_barcodes, secondary_barcodes],
            output_directory = output_directory,
            sample_key = sample_key,
            sample_list = sample_list,
            read_count = self.file_extractor.get_read_count(),
            parse_type = parse_type)
        
        qseq_file_parser.write()    