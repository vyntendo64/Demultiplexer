#!/usr/bin/env python3

from MultipleIterator import MultipleSequencingFileIterator

class QseqFileParser:
    def __init__(self, files):
        self.files = files
        self.output_dict = {}

    def run(self):
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
            barcode_indexes = self.duplicates(self.action, 'barcode')
            # get position of read files
            read_indexes = self.duplicates(self.action, 'read')
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
                            output_object.write(self.qseq_fastq_conversion(line[read_indexes[out_count]]))
                    else:
                        # if barcode sequence not in dictionary write to unmatched
                        self.unmatched_read += 1
                        sample = 'unmatched'
                        out = self.output_dict[sample]
                        for out_count, output_object in enumerate(out):
                            output_object.write(self.qseq_fastq_conversion(line[read_indexes[out_count]]))
        # close all output objects
        for sample in self.output_dict.values():
            for out_object in sample:
                out_object.close()

    def duplicates(self, lst, item):
        """Python index lookup in a list returns the first index by default, this function returns all indexes of an item
        in a list. Loops through list and returns count if item = item in list.
           -----------------------------------------------------
           lst=List
           item=object in list
           returns; list with indices of object"""
        return [i for i, x in enumerate(lst) if x == item]


    def qseq_fastq_conversion(self, qseq_list):
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