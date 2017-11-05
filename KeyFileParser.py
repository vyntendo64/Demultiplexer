#!/usr/bin/env python3

class KeyFileParser:
    """Opens Illumina qseq directory and processes qseq files, outputs samples fastq files"""

    def __init__(self, path = 'path', parse_type = 'single'):
        self.path = path
        self.parse_type = parse_type
        self.sample_list = []

    def set_parse_type(self, parse_type = 'single'):
        self.parse_type = parse_type

    def get_parse_type(self):
        return self.parse_type

    def get_sample_list(self):
        """Takes sample label file and processes barcode sample IDs.  Note this function assumes 'barcode1 \t barcode 2
        \t sample_name \n' .  If only one barcode is used then the sample file should be formatted as
        'barcode1 \t sample_name \n' and no barcode_2 file should be supplied. Function will fail if sample key is not
        in this format.
        -----------------------------------------------------
        opens self.sample_key; parses files and hashes to sample name based on unique barcode ID
        returns self.sample_key; a dictionary hashing to sample name"""
        # initialize dict
        # loop through text file

        sample_list = []
        for line in open(self.path):
            # replace new line indicator
            line = line.replace('\n', '').split('\t')
            # if barcode_2 supplied sample id is a combination of 2 barcodes
            try:
                if self.parse_type == 'dual':
                    sample_list.append({line[1]: [int(line[0]), int(line[1])]})
                else: 
                    sample_list.append({line[1]: [int(line[0])]})
            except ValueError:
                print('Sample Key not in proper format\nPlease format file Barcode tab SampleName')
                sys.exit()

        return sample_list


