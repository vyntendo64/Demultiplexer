#!/usr/bin/env python3

class KeyFileParser:
    """Opens Illumina qseq directory and processes qseq files, outputs samples fastq files"""

    def __init__(self, path = 'path'):
        self.path = path
        self.sample_list = []
        self.sample_key = {}

    def get_sample_list(self):
        return self.sample_list

    def get_sample_key(self):
        return self.sample_key

    def set_labels(self, to = 'single'):
        if to == 'dual':
            self.set_dual_labels()

        else:
            self.set_single_labels()

    def set_dual_labels(self):
        """Takes sample label file and processes barcode sample IDs.  Note this function assumes 'barcode1 \t barcode 2
        \t sample_name \n' .  If only one barcode is used then the sample file should be formatted as
        'barcode1 \t sample_name \n' and no barcode_2 file should be supplied. Function will fail if sample key is not
        in this format.
        -----------------------------------------------------
        opens self.sample_key; parses files and hashes to sample name based on unique barcode ID
        returns self.sample_key; a dictionary hashing to sample name"""
        # initialize dict
        
        sample_dict = {}
        # loop through text file
        for line in open(self.path):
            # replace new line indicator
            line_replace = line.replace('\n', '')
            # split on tabs
            line_split = line_replace.split('\t')
            # if barcode_2 supplied sample id is a combination of 2 barcodes
            try:
                sample_dict['key' + str(int(line_split[0])) + 'key' + str(int(line_split[1]))] = line_split[2]
            except ValueError:
                print('Sample Key not in proper format\nPlease format file Barcode1 tab Barcode2 tab SampleName')
                sys.exit()

            self.sample_list.append(line_split[2])
        self.sample_key = sample_dict

    def set_single_labels(self):
        """Takes sample label file and processes barcode sample IDs.  Note this function assumes 'barcode1 \t barcode 2
        \t sample_name \n' .  If only one barcode is used then the sample file should be formatted as
        'barcode1 \t sample_name \n' and no barcode_2 file should be supplied. Function will fail if sample key is not
        in this format.
        -----------------------------------------------------
        opens self.sample_key; parses files and hashes to sample name based on unique barcode ID
        returns self.sample_key; a dictionary hashing to sample name"""
        # initialize dict
        # loop through text file

        sample_dict = {}
        for line in open(self.path):
            # replace new line indicator
            line_replace = line.replace('\n', '')
            # split on tabs
            line_split = line_replace.split('\t')
            # if barcode_2 supplied sample id is a combination of 2 barcodes
            try:
                sample_dict['key' + str(int(line_split[0]))] = line_split[1]
            except ValueError:
                print('Sample Key not in proper format\nPlease format file Barcode tab SampleName')
                sys.exit()
            self.sample_list.append(line_split[1])
        self.sample_key = sample_dict