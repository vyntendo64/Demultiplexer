#!/usr/bin/env python3

from BarcodeFileParser import BarcodeFileParser
from DirectoryFileParser import DirectoryFileParser
from KeyFileParser import KeyFileParser
from QseqFileParser import QseqFileParser

import sys

class Demuliplex:
    """Opens Illumina qseq directory and processes qseq files, outputs samples fastq files"""

    def __init__(self, 
        *args, 
        directory = 'path',
        output_directory = 'path', 
        sample_key = 'path', 
        mismatch = 1, 
        action = '', 
        barcode_1 = None,
        barcode_2 = None, 
        gnu_zipped = False):

        #HUH
        # store file description
        self.file_description = []

        for arg in args:
            self.file_description.append(arg.split('^'))

        # check if all input files have labels
        if len(action) != len(self.file_description):
            #print('# of input files not equal to the number of input file labels')
            sys.exit()

        # set input variables
        self.directory = directory
        self.output_directory = output_directory
        self.sample_key = sample_key
        self.mismatch = mismatch
        self.action = action
        self.barcode_1 = barcode_1
        self.barcode_2 = barcode_2
        self.gnu_zipped = gnu_zipped

        self.setUp()
        
        directory_parser = DirectoryFileParser(directory = directory, 
            file_description = self.file_description)

        self.files = directory_parser.get_files()


    def setUp(self):
        self.files = []
        self.barcode_count = None
        self.read_count = None
        self.reads = 0
        self.reads_pass_filter = 0
        self.unmatched_read = 0
        self.indexed_reads = 0
        self.sample_list = []
        self.output_dict = {}

    def run(self):
        # print('sself filess')
        # print(self.files)
        self.process_barcodes()
        self.set_action()

        key_file_parser = KeyFileParser(path = self.sample_key)

        if self.barcode_2:
            key_file_parser.set_combination_labels()
        else:
            key_file_parser.set_single_labels()

        self.sample_list = key_file_parser.get_sample_list();
        self.sample_key = key_file_parser.get_sample_key();

        qseq_file_parser = QseqFileParser(files = self.files,
            directory = self.directory,
            gnu_zipped = self.gnu_zipped)
        
        qseq_file_parser.run()


    def process_barcodes(self):
        """If barcode file supplied process files and store values in dictionary
        -----------------------------------------------------
        opens self.barcode*, a path to text file with a new barcode on each line
        returns self.barcode*, a dictionary hashing barcodes to an Illumina ID"""
        # todo broke tests by changing hash collision detection

        barcode_parser = BarcodeFileParser(self.barcode_1)
        self.barcode_1 = barcode_parser.get_barcodes()

        barcode_parser = BarcodeFileParser(self.barcode_2)
        self.barcode_2 = barcode_parser.get_barcodes()

    def set_action(self):
        """Parses string describing input files, the action should be formatted as r for read and b for barcode,
        so a string 'rbbr' would describe a four qseq file input with read barcode barcode read.  Used to properly
        parse the qseq files
        -----------------------------------------------------
        self.action: string describing input files
        self.barcode_count: returns count of barcode in file_lablel as a downstream control"""

        label_list = []
        for character in self.action:
            if character.lower() == 'r':
                label_list.append('read')
            elif character.lower() == 'b':
                label_list.append('barcode')

        self.action = label_list
        self.barcode_count = label_list.count('barcode')
        self.read_count = label_list.count('read')

    def get_sample_labels(self):
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
        for line in open(self.sample_key):
            # replace new line indicator
            line_replace = line.replace('\n', '')
            # split on tabs
            line_split = line_replace.split('\t')
            # if barcode_2 supplied sample id is a combination of 2 barcodes
            if self.barcode_2:
                # int(), to check if file in is proper format
                try:
                    sample_dict['key' + str(int(line_split[0])) + 'key' + str(int(line_split[1]))] = line_split[2]
                except ValueError:
                    # print('Sample Key not in proper format\nPlease format file Barcode1 tab Barcode2 tab SampleName')
                    sys.exit()
                self.sample_list.append(line_split[2])
            # else id is only one barcode
            else:
                try:
                    sample_dict['key' + str(int(line_split[0]))] = line_split[1]
                except ValueError:
                    # print('Sample Key not in proper format\nPlease format file Barcode tab SampleName')
                    sys.exit()
                self.sample_list.append(line_split[1])
        self.sample_key = sample_dict

    

    def output_objects(self):
        """Initialized objects to output reads in fastq format, will generate a file for every 'read' labeled file in
        the action plus a file for unmatched reads
        -----------------------------------------------------
        output_directory = path to write files; folder must already exist
        self.sample_list: list of input samples
        self.read_count: number of read files labeled in file label
        returns self.output_dict; hashes to output object based on sample name"""
        # initialize output objects for all samples
        # print(self.sample_list)
        # print('read count')
        # print(self.read_count)
        for sample in self.sample_list:
            object_list = []
            for count in range(self.read_count):
                object_list.append(open(self.output_directory + sample + '_' + str(count + 1) + '.fastq', 'w'))
            self.output_dict[sample] = object_list
        object_list = []
        # initialize output objects for unmatched reads
        for count in range(self.read_count):
            object_list.append(open(self.output_directory + 'unmatched' + '_' + str(count + 1) + '.fastq', 'w'))

        self.output_dict['unmatched'] = object_list

    def iterate_through_qseq(self):
        """Iterate through groups of ID'd qseq files and output demultiplexed fastq files
        -----------------------------------------------------
        self.files = sorted list of input files;
        self.barcode*: dictionary of barcodes
        self.action: list of barcode/ read file positions
        self.output_dict: dict of output objects, 1 per read qseq files
        returns;
        self.read = number of reads processes across all of the groups of files
        self.reads_pass_filter = number of reads that pass the Illumina filter
        self.indexed_reads = reads matched to sample index
        self.unmatched_reads = number of unmatched reads
        """
        # transpose iterator list
        self.files = list(map(list, zip(*self.files)))
        # loop through lists of files
        for files in self.files:
            # initialize iterator object for sorted group of files
            iterator = MultipleSequencingFileIterator(*files, directory=self.directory, gnu_zipped=self.gnu_zipped)
            # get position of barcode files
            barcode_indexes = duplicates(self.action, 'barcode')
            # get position of read files
            read_indexes = duplicates(self.action, 'read')
            # set barcode list, for looping
            barcode_list = [self.barcode_1, self.barcode_2]
            # loop through grouped files
            for count, line in enumerate(iterator.iterator_zip()):
                self.reads += 1
                # set string with Illumina quality control information
                combined_filter = ''.join([qual[-1] for qual in line])
                # initialize empty sample_id value
                sample_id = ''
                # if all reads don't pass filter don't consider
                if '0' not in combined_filter:
                    self.reads_pass_filter += 1
                    # loop through barcode_indexes, get sample key
                    for index_count, index in enumerate(barcode_indexes):
                        try:
                            # get sequence location in qseq file
                            key = barcode_list[index_count][line[barcode_indexes[index_count]][8]]
                        except KeyError:
                            # if barcode sequence not in barcode dictionary set key to 'x'
                            key = 'x'
                        sample_id = '{0}key{1}'.format(sample_id, str(key))
                    # if barcode matches with key proceed
                    if 'x' not in sample_id:
                        try:
                            # look up sample, if matched set sample name
                            sample = self.sample_key[sample_id]
                            self.indexed_reads += 1
                        except KeyError:
                            # if sample unmatched write to unmatched reads
                            self.unmatched_read += 1
                            sample = 'unmatched'
                        # retrieve list of output objects
                        out = self.output_dict[sample]
                        # write line to file
                        for out_count, output_object in enumerate(out):
                            # convert qseq line to fatq format
                            output_object.write(qseq_fastq_conversion(line[read_indexes[out_count]]))
                    else:
                        # if barcode sequence not in dictionary write to unmatched
                        self.unmatched_read += 1
                        sample = 'unmatched'
                        out = self.output_dict[sample]
                        for out_count, output_object in enumerate(out):
                            output_object.write(qseq_fastq_conversion(line[read_indexes[out_count]]))
        # close all output objects
        for sample in self.output_dict.values():
            for out_object in sample:
                out_object.close()
