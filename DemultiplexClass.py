#!/usr/bin/env python3

from MultipleIterator import MultipleSequencingFileIterator
from os import listdir
import sys


def reverse_complement(string):
    """Simple reverse complement function used to initialize a barcode dictionary, (original sequence and reverse
    complement both link to same hash in dictionary).
    -----------------------------------------------------
    string='string': string must be composed of BP ATGC
    returns; reverse complement string"""
    # reverse string
    reverse_string = string[::-1]
    # complementary bp lookup dictionary
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    all_bases = list(reverse_string)
    complement_list = []
    # iterate over string list
    for i in all_bases:
        complement_list.append(complement[i])
    # return joined string
    return ''.join(complement_list)


def duplicates(lst, item):
    """Python index lookup in a list returns the first index by default, this function returns all indexes of an item
    in a list. Loops through list and returns count if item = item in list.
       -----------------------------------------------------
       lst=List
       item=object in list
       returns; list with indices of object"""
    return [i for i, x in enumerate(lst) if x == item]


def qseq_fastq_conversion(qseq_list):
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


def barcode_mismatch(barcode='string', mismatch_list=('A', 'T', 'G', 'C', '.'), number_mismatches=2):
    """Takes a barcode input and outputs a list containing the original barcode and mistmatched barcodes. Will return
    every possible mismatch based on the mismatch list input.  Slow implementation of designating mismatches, but
    the script is only used to initialize a hash table. This method is much faster then calculating the hamming distance
    for every barcode (just pre-computing mismatched barcodes).
    -----------------------------------------------------
    barcode;'string' typical Illumina barcodes are 8 bp
    mismatch_list;list of mismatch characters to be inserted at every barcode position
    number_mismatches;int number of possible mismatches compared to original barcode to compute, anything above 2
    mismatches decreases demultiplexing performance
    returns; a list of barcodes ('strings')"""
    # initialize list with with input barcode
    all_barcodes = [barcode]
    # calculate first set of mismatches
    if number_mismatches != 0:
        # loop over possible mismatches
        for possible_mismatch in mismatch_list:
            # loop over every position in barcode
            for character in range(len(barcode)):
                # barcode to list
                barcode_list = list(barcode)
                # set character in to mismatch
                barcode_list[character] = possible_mismatch
                # if mismatched barcode not already in list then add the barcode to the list
                if ''.join(barcode_list) not in all_barcodes:
                    all_barcodes.append(''.join(barcode_list))
        if number_mismatches > 1:
            # initialize additional loops based on mismatch number
            for count in range(number_mismatches - 1):
                # repeat loop
                mismatched_barcodes = list(all_barcodes)
                for mismatched_barcode in mismatched_barcodes:
                    for possible_mismatch in mismatch_list:
                        for character in range(len(mismatched_barcode)):
                            barcode_list = list(mismatched_barcode)
                            barcode_list[character] = possible_mismatch
                            if ''.join(barcode_list) not in all_barcodes:
                                all_barcodes.append(''.join(barcode_list))
    return list(set(all_barcodes))


class Demuliplex:
    """Opens Illumina qseq directory and processes qseq files, outputs samples fastq files"""

    def __init__(self, *args, directory='path', sample_key='path', mismatch=1, file_label='', barcode_1=None,
                 barcode_2=None, gnu_zipped=False):
        # store file description
        self.file_description = []
        for arg in args:
            self.file_description.append(arg.split('^'))
        # check if all input files have labels
        if len(file_label) != len(self.file_description):
            print('# of input files not equal to the number of input file labels')
            sys.exit()
        # set input variables
        self.directory = directory
        self.mismatch = mismatch
        self.file_label = file_label
        self.barcode_1 = barcode_1
        self.barcode_2 = barcode_2
        self.sample_key = sample_key
        self.file_list = []
        self.output_dict = {}
        self.barcode_count = None
        self.read_count = None
        self.reads = 0
        self.reads_pass_filter = 0
        self.unmatched_read = 0
        self.indexed_reads = 0
        self.sample_list = []
        self.gnu_zipped = gnu_zipped

    def process_barcodes(self):
        """If barcode file supplied process files and store values in dictionary
        -----------------------------------------------------
        opens self.barcode*, a path to text file with a new barcode on each line
        returns self.barcode*, a dictionary hashing barcodes to an Illumina ID"""
        # todo broke tests by changing hash collision detection
        if self.barcode_1:
            # initialize dict
            barcode_dict = {}
            # loop through barcode text file
            for count, line in enumerate(open(self.barcode_1)):
                # replace new line indicator
                barcode = line.replace('\n', '')
                # look up reverse complement
                reverse_barcode = reverse_complement(barcode)
                # get list containing barcode and all mismatches
                barcode_list = barcode_mismatch(barcode, number_mismatches=self.mismatch)
                # iterate over list and initialize a dictionary key value pair for every barcode
                for possible_barcode in barcode_list:
                    if possible_barcode not in list(barcode_dict.keys()):
                        barcode_dict[possible_barcode] = count + 1
                    else:
                        if count + 1 != barcode_dict[possible_barcode]:
                            print('Hash Collision, reduce the number of allowed mismatches')
                            sys.exit()
                # repeat for barcode reverse complement
                reverse_barcode_list = barcode_mismatch(reverse_barcode, number_mismatches=self.mismatch)
                for possible_reverse_barcode in reverse_barcode_list:
                    if possible_reverse_barcode not in list(barcode_dict.keys()):
                        barcode_dict[possible_reverse_barcode] = count + 1
                    else:
                        if count + 1 != barcode_dict[possible_reverse_barcode]:
                            print('Hash Collision, reduce the number of allowed mismatches')
                            sys.exit()
            self.barcode_1 = barcode_dict
        # if barcode_2 file is given, repeat process
        if self.barcode_2:
            barcode_dict = {}
            for count, line in enumerate(open(self.barcode_2)):
                # replace new line indicator
                barcode = line.replace('\n', '')
                # look up reverse complement
                reverse_barcode = reverse_complement(barcode)
                barcode_list = barcode_mismatch(barcode, number_mismatches=self.mismatch)
                for possible_barcode in barcode_list:
                    if possible_barcode not in list(barcode_dict.keys()):
                        barcode_dict[possible_barcode] = count + 1
                    else:
                        if count + 1 != barcode_dict[possible_barcode]:
                            print('Hash Collision, reduce the number of allowed mismatches')
                            sys.exit()
                reverse_barcode_list = barcode_mismatch(reverse_barcode, number_mismatches=self.mismatch)
                for possible_reverse_barcode in reverse_barcode_list:
                    if possible_reverse_barcode not in list(barcode_dict.keys()):
                        barcode_dict[possible_reverse_barcode] = count + 1
                    else:
                        if count + 1 != barcode_dict[possible_reverse_barcode]:
                            print('Hash Collision, reduce the number of allowed mismatches')
                            sys.exit()
            self.barcode_2 = barcode_dict

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
                    print('Sample Key not in proper format\nPlease format file Barcode1 tab Barcode2 tab SampleName')
                    sys.exit()
                self.sample_list.append(line_split[2])
            # else id is only one barcode
            else:
                try:
                    sample_dict['key' + str(int(line_split[0]))] = line_split[1]
                except ValueError:
                    print('Sample Key not in proper format\nPlease format file Barcode tab SampleName')
                    sys.exit()
                self.sample_list.append(line_split[1])
        self.sample_key = sample_dict

    def get_directory_lists(self):
        """Link to directory, pull  list of files, combine files with same unique id. Function only works with files
        that have an integer as a unique id
        -----------------------------------------------------
        self.directory: path to illumina sequencing lane directory
        self.file_description: list of lists containing sequencing file prefix and suffix, sequencing file designated
        by prefix*suffix
        returns; sorted list of relevant files names in directory"""
        file_list = listdir(self.directory)
        # initialize list to hold sample names (ie. coupled file IDs)
        sample_names = [[] for _ in range(len(self.file_description))]
        # key to sort files bases on proper ID
        sorting_key = [[] for _ in range(len(self.file_description))]
        for count, file_title in enumerate(self.file_description):
            for file_name in file_list:
                # order files based on read type, ie all read 1 go together etc.
                if file_title[0] in file_name:
                    sample_names[count].append(file_name)
                    # store unique file ID
                    try:
                        sort_id = int((file_name.replace(file_title[0], '')).replace(file_title[1], ''))
                    except ValueError:
                        print('File ID must be integer')
                        sys.exit()
                    # store sort id in list for file type
                    sorting_key[count].append(sort_id)
        for count, seq_file_list in enumerate(sample_names):
            # sort list on unique ID and append sorted list of file, one list per file prefix
            self.file_list.append([x for x, y in sorted(zip(sample_names[count], sorting_key[count]))])

    def process_file_label(self):
        """Parses string describing input files, the file_label should be formatted as r for read and b for barcode,
        so a string 'rbbr' would describe a four qseq file input with read barcode barcode read.  Used to properly
        parse the qseq files
        -----------------------------------------------------
        self.file_label: string describing input files
        self.barcode_count: returns count of barcode in file_lablel as a downstream control"""
        label_list = []
        for character in self.file_label:
            if character.lower() == 'r':
                label_list.append('read')
            elif character.lower() == 'b':
                label_list.append('barcode')
        self.file_label = label_list
        self.barcode_count = label_list.count('barcode')
        self.read_count = label_list.count('read')

    def output_objects(self, output_directory='path'):
        """Initialized objects to output reads in fastq format, will generate a file for every 'read' labeled file in
        the file_label plus a file for unmatched reads
        -----------------------------------------------------
        output_directory = path to write files; folder must already exist
        self.sample_list: list of input samples
        self.read_count: number of read files labeled in file label
        returns self.output_dict; hashes to output object based on sample name"""
        # initialize output objects for all samples
        for sample in self.sample_list:
            object_list = []
            for count in range(self.read_count):
                object_list.append(open(output_directory + sample + '_' +
                                        str(count + 1) + '.fastq', 'w'))
            self.output_dict[sample] = object_list
        object_list = []
        # initialize output objects for unmatched reads
        for count in range(self.read_count):
            object_list.append(open(output_directory + 'unmatched' + '_' +
                                    str(count + 1) + '.fastq', 'w'))
        self.output_dict['unmatched'] = object_list

    def iterate_through_qseq(self):
        """Iterate through groups of ID'd qseq files and output demultiplexed fastq files
        -----------------------------------------------------
        self.file_list = sorted list of input files;
        self.barcode*: dictionary of barcodes
        self.file_label: list of barcode/ read file positions
        self.output_dict: dict of output objects, 1 per read qseq files
        returns;
        self.read = number of reads processes across all of the groups of files
        self.reads_pass_filter = number of reads that pass the Illumina filter
        self.indexed_reads = reads matched to sample index
        self.unmatched_reads = number of unmatched reads
        """
        # transpose iterator list
        self.file_list = list(map(list, zip(*self.file_list)))
        # loop through lists of files
        for files in self.file_list:
            # initialize iterator object for sorted group of files
            iterator = MultipleSequencingFileIterator(*files, directory=self.directory, gnu_zipped=self.gnu_zipped)
            # get position of barcode files
            barcode_indexes = duplicates(self.file_label, 'barcode')
            # get position of read files
            read_indexes = duplicates(self.file_label, 'read')
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
