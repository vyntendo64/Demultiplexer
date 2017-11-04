#!/usr/bin/env python3
import os
from MultipleIterator import MultipleSequencingFileIterator

class QseqFileParser:
    def __init__(self, 
        files = [], 
        barcode_list = [], 
        output_directory = os.getcwd(), 
        sample_key = {},
        sample_list = [],
        read_count = 1,
        gnu_zipped = False): 
        
        self.files = files
        self.barcode_list = barcode_list
        self.output_directory = output_directory
        self.sample_key = sample_key
        self.sample_list = sample_list
        self.read_count = read_count
        self.gnu_zipped = gnu_zipped

    def write(self):
        self.reads = 0
        self.reads_pass_filter = 0
        self.unmatched_reads = 0
        self.indexed_reads = 0
        self.UNMATCHED_SAMPLE = 'unmatched'

        output_dict = self.get_output_dict()
        print('output_dict')
        print(output_dict)
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
        print('running')
        print(self.files)
        print('zip files')
        print(zip(*self.files));
        print('zip files')
        print(map(list, zip(*self.files)));
        print('zip files')
        print(list(map(list, zip(*self.files))));
        unizipped_files = list(map(list, zip(*self.files)))
        print(unizipped_files)
        # loop through lists of files
        for file in self.files:
            print('file')
            print(file)
            # print('files')
            # print(files)
            # initialize iterator object for sorted group of files
            iterator = MultipleSequencingFileIterator(files = files, 
                gnu_zipped = self.gnu_zipped)
            print('iterator')
            print(iterator)

            # get position of barcode files
            barcode_indexes = self.get_duplicates_count(self.action, 'barcode')
            # print(barcode_indexes)
            # print('barcode_indexes')
            # get position of read files
            read_indexes = self.get_duplicates_count(self.action, 'read')
            # set barcode list, for looping
            # loop through grouped files
            print('zip')
            print(enumerate(iterator.iterator_zip()))
            for count, line in enumerate(iterator.iterator_zip()):
                print('count')
                print(count)
                print('line')
                print(line)

                sample_id = ''
                self.increment_reads()
                # set string with Illumina quality control information
                combined_filter = self.get_combined_filter(line = line)
                print('combined filter')
                print(combined_filter)
                # initialize empty sample_id value
                # print(combined_filter)
                # print('combined_filter')
                # if all reads don't pass filter don't consider
                if '0' not in combined_filter:
                    self.increment_reads_pass_filter()
                    # loop through barcode_indexes, get sample key
                    sample_id = self.get_sample_id(barcode_indexes = barcode_indexes, line = line)
                    print('sample_id')
                    print(sample_id)
                    # print(sample_id)    
                    # print('sample_id')                
                    # if barcode matches with key proceed
                    if 'x' not in sample_id:
                        try:
                            # look up sample, if matched set sample name
                            sample = self.sample_key[sample_id]
                            self.increment_indexed_reads()
                        except KeyError:
                            # if sample unmatched write to unmatched reads
                            self.increment_unmatched_reads()
                            sample = 'unmatched'
                        # retrieve list of output objects
                        out = self.output_dict[sample]
                        # write line to file
                        for out_count, output_object in enumerate(out):
                            # convert qseq line to fatq format
                            print('out_count')
                            print(out_count)
                            print('read_indexes')
                            print(read_indexes)
                            print('output object')
                            print(output_object)

                            output_object.write(self.qseq_fastq_conversion(line[read_indexes[out_count]]))
                    else:
                        # if barcode sequence not in dictionary write to unmatched
                        self.increment_unmatched_reads()
                        sample = 'unmatched'
                        out = self.output_dict[sample]
                        for out_count, output_object in enumerate(out):
                            output_object.write(self.qseq_fastq_conversion(line[read_indexes[out_count]]))
        # close all output objects
        for sample in self.output_dict.values():
            for out_object in sample:
                out_object.close()

    def get_output_dict(self):
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
        output_dict = {}
        self.sample_list.append('unmatched')
        for action in self.sample_list:
            object_list = []
            for count in range(self.read_count):
                object_list.append(open(self.output_directory + '/' + action + '_' + str(count + 1) + '.fastq', 'w'))
            output_dict[action] = object_list
        return output_dict

    def get_sample_id(self, barcode_indexes, line):
        sample_id = ''
        for index_count, index in enumerate(barcode_indexes):
            try:
                # get sequence location in qseq file
                key = self.barcode_list[index_count][line[barcode_indexes[index_count]][8]]
            except KeyError:
                # if barcode sequence not in barcode dictionary set key to 'x'
                key = 'x'
            sample_id = '{0}key{1}'.format(sample_id, str(key))

        return sample_id

    def increment_reads(self): 
        self.reads += 1

    def increment_reads_pass_filter(self):
        self.reads_pass_filter += 1

    def increment_unmatched_reads(self):
        self.unmatched_reads += 1

    def increment_indexed_reads(self):
        self.indexed_reads += 1

    def get_combined_filter(self, line):
        return ''.join([qual[-1] for qual in line])

    def get_duplicates_count(self, lst, item):
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

    