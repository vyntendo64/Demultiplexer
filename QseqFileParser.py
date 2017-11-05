#!/usr/bin/env python3
import os
import sys
from MultipleIterator import MultipleSequencingFileIterator

class QseqFileParser:
    def __init__(self, 
        files = [], 
        barcode_list = [], 
        output_directory = os.getcwd(), 
        sample_list = [],
        read_count = 1): 
        
        self.files = files
        self.barcode_list = barcode_list
        self.output_directory = output_directory
        self.sample_list = sample_list
        self.read_count = read_count
        
        self.BARCODE_LOCATION = 8

    def run(self, 
        files = [], 
        gnu_zipped = False):

        self.reads = 0
        self.reads_pass = 0
        self.unmatched_reads = 0
        self.indexed_reads = 0

        self.iterator = MultipleSequencingFileIterator(files = files, gnu_zipped = gnu_zipped)

        self.barcode_indexes = self.iterator.get_barcode_indexes()
        self.read_indexes = self.iterator.get_read_indexes()

        output_dict = self.get_output_dict()

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

        for line_no, line in enumerate(list(self.iterator.get_next_row())):
            self.increment_reads()

            if self.line_is_safe(line):
                self.increment_reads_pass()
                sample_ids = self.get_sample_ids(line = line)
                sample = self.get_sample(sample_ids)

                if self.sample_is_ok(sample):
                    self.increment_indexed_reads()
                    self.write(output_dict = output_dict, sample = sample, line = line)
                    # retrieve list of output objects
                else:
                    # if barcode sequence not in dictionary write to unmatched
                    self.increment_unmatched_reads()
                    self.write(output_dict = output_dict, sample = 'unmatched', line = line)

        # close all output objects
        self.close(output_dict)

    def get_output_dict(self):
        """Initialized objects to output reads in fastq format, will generate a file for every 'read' labeled file in
        the action plus a file for unmatched reads
        -----------------------------------------------------
        output_directory = path to write files; folder must already exist
        self.sample_list: list of input samples
        self.read_count: number of read files labeled in file label
        returns self.output_dict; hashes to output object based on sample name"""
        # initialize output objects for all samples
        output_dict = {}
        self.sample_list.append({'unmatched': []})
        for sample in self.sample_list:
            for sample_no in sample:
                object_list = []
                for count in range(self.read_count):
                    object_list.append(open(self.output_directory + '/' + sample_no + '_' + str(count + 1) + '.fastq', 'w'))
                output_dict[sample_no] = object_list
        return output_dict

    def line_is_safe(self, line):
        for qual in line:
            if qual[-1] != "1":
                return False
        return True    

    def sample_ids_pass(self, sample_ids):
        samples_status = True

        for sample_id in sample_ids:
            if sample_id == 'x':
                samples_status = False

        return samples_status

    def get_sample_ids(self, line):
        sample_ids = []

        for barcode_file_no, barcode_index in enumerate(self.barcode_indexes):
            sample_id = self.get_sample_id(line = line, barcode_file_no = barcode_file_no)
            sample_ids.append(sample_id)
        return sample_ids

    def get_sample_id(self, line, barcode_file_no):
        try:
            # get sequence location in qseq file
            key = self.barcode_list[barcode_file_no][line[self.barcode_indexes[barcode_file_no]][self.BARCODE_LOCATION]]
        except KeyError:
            # if barcode sequence not in barcode dictionary set key to 'x'
            key = 'x'
        return key

    def get_sample(self, sample_ids):
        return_sample = 'unmatched'
        for sample_item in self.sample_list:
            for sample in sample_item:
                if sample_item[sample] == sample_ids:
                    return_sample = sample
                    break
        return return_sample

    def sample_is_ok(self, sample):
        if(sample == 'unmatched'):
            return False
        return True

    def write(self, output_dict, sample, line):
        out = output_dict[sample]
        # write line to file
        for out_count, output_object in enumerate(out):
            # convert qseq line to fatq format
            output_object.write(self.format(line[self.read_indexes[out_count]]))

    def format(self, qseq_list):
        """Convert a .qseq file to .fastq before output
           -----------------------------------------------------
           qseq_list: '\t' split qseq string
           returns; an output string with four new line indicators in fastq format"""
        # pull fastq header information
        fastq_id = '@%s:%s:%s:%s:%s#%s/%s' % (qseq_list[0], qseq_list[2], qseq_list[3], qseq_list[4],
                                              qseq_list[5], qseq_list[6], qseq_list[7])
        # mask any missing base calls
        seq = qseq_list[8].replace('.', 'N')
        # line 3 is used by some tools to store additional information, so the output will be static
        line_3 = '+'
        # The quality output is a single value in qseq list
        quality = qseq_list[9]

        fastq_out = fastq_id + '\n' + seq + '\n' + line_3 + '\n' + quality + '\n'

        return fastq_out

    def close(self, output_dict):
        for sample in output_dict.values():
            for out_object in sample:
                out_object.close()

    def increment_reads(self): 
        self.reads += 1

    def increment_reads_pass(self):
        self.reads_pass += 1

    def increment_unmatched_reads(self):
        self.unmatched_reads += 1

    def increment_indexed_reads(self):
        self.indexed_reads += 1