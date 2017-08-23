#! /usr/env python3

from MultipleIterator import MultipleSequencingFileIterator
from os import listdir


def hamming_distance(s1, s2):
    """Used to calculate mismatches between barcode and sequence call"""
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def reverse_complement(string):
    reverse_string = string[::-1]
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    all_bases = list(reverse_string)
    complement_list = []
    for i in all_bases:
        complement_list.append(complement[i])
    return ''.join(complement_list)


class Demuliplex:

    def __init__(self, *args,  directory='path', sample_key='path', mismatch=1, file_label='rbbr', barcode_1=None, barcode_2=None):
        self.file_description = []
        for arg in args:
            self.file_description.append(arg.split('*'))
        if len(file_label) != len(self.file_description):
            print('something is wrong')
        self.directory = directory
        self.mismatch = mismatch
        self.file_label = file_label
        self.barcode_1 = barcode_1
        self.barcode_2 = barcode_2
        self.sample_key = sample_key
        self.file_list = []

    def process_barcodes(self):
        if self.barcode_1:
            barcode_dict = {}
            for count, line in enumerate(open(self.barcode_1)):
                barcode = line.replace('\n', '')
                reverse_barcode = reverse_complement(barcode)
                barcode_dict[count+1] = [barcode, reverse_barcode]
            self.barcode_1 = barcode_dict
        if self.barcode_2:
            barcode_dict = {}
            for count, line in enumerate(open(self.barcode_2)):
                barcode = line.replace('\n', '')
                reverse_barcode = reverse_complement(barcode)
                barcode_dict[count+1] = [barcode, reverse_barcode]
            self.barcode_2 = barcode_dict

    def get_sample_labels(self):
        """Barcode1, barcode2, sample_name"""
        sample_dict = {}
        for line in open(self.sample_key):
            line_replace = line.replace('\n', '')
            line_split = line_replace.split('\t')
            sample_dict[line_split[2]] = [int(line[0]), int(line[1])]
        self.sample_key = sample_dict

    def get_directory_lists(self):
        file_list = listdir(self.directory)
        sample_names = [[] for x in range(len(self.file_description))]
        sorting_key = [[] for y in range(len(self.file_description))]
        for count, file_title in enumerate(self.file_description):
            for file_name in file_list:
                if file_title[0] in file_name:
                    sample_names[count].append(file_name)
                    sort_id = int((file_name.replace(file_title[0], '')).replace(file_title[1], ''))
                    sorting_key[count].append(sort_id)
        for count, seq_file_list in enumerate(sample_names):
            self.file_list.append([x for x, y in sorted(zip(sample_names[count], sorting_key[count]))])
        print(self.file_list)


    #def output_fastq:

    #def iterate_through_fastq:


test = Demuliplex('s_1_4_*_qseq.txt', 's_1_2_*_qseq.txt', directory='/home/colin/Dropbox/Pelligrini_Lab')
test.get_directory_lists()