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


def duplicates(lst, item):
    return [i for i, x in enumerate(lst) if x == item]


def qseq_fastq_conversion(qseq_list):
    fastq_id = '@%s:%s:%s:%s:%s#%s/%s' % (qseq_list[0], qseq_list[2], qseq_list[3], qseq_list[4],
                                          qseq_list[5], qseq_list[6], qseq_list[7])
    seq = qseq_list[8].replace('.', 'N')
    line_3 = '+'
    quality = qseq_list[9]
    fastq_out = fastq_id + '\n' + seq + '\n' + line_3 + '\n' + quality + '\n'
    return fastq_out


class Demuliplex:

    def __init__(self, *args, directory='path', sample_key='path', mismatch=1, file_label='rbbr', barcode_1=None,
                 barcode_2=None):
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
        self.output_dict = {}
        self.barcode_count = None
        self.reads = 0
        self.reads_pass_filter = 0
        self.unmatched_read = 0
        self.sample_list = []

    def process_barcodes(self):
        if self.barcode_1:
            barcode_dict = {}
            for count, line in enumerate(open(self.barcode_1)):
                barcode = line.replace('\n', '')
                reverse_barcode = reverse_complement(barcode)
                barcode_dict[barcode] = count + 1
                barcode_dict[reverse_barcode] = count + 1
            self.barcode_1 = barcode_dict
        if self.barcode_2:
            barcode_dict = {}
            for count, line in enumerate(open(self.barcode_2)):
                barcode = line.replace('\n', '')
                reverse_barcode = reverse_complement(barcode)
                barcode_dict[barcode] = count + 1
                barcode_dict[reverse_barcode] = count + 1
            self.barcode_2 = barcode_dict

    def get_sample_labels(self):
        """Barcode1, barcode2, sample_name"""
        sample_dict = {}
        for line in open(self.sample_key):
            line_replace = line.replace('\n', '')
            line_split = line_replace.split('\t')
            if self.barcode_2:
                sample_dict['key' + str(int(line_split[0])) + ' key' + str(int(line_split[1]))] = line_split[2]
                self.sample_list.append(line_split[2])
            else:
                sample_dict[line_split[1]] = int(line_split[0])
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

    def process_file_label(self):
        label_list = []
        for character in self.file_label:
            if character.lower() == 'r':
                label_list.append('read')
            elif character.lower() == 'b':
                label_list.append('barcode')
        self.file_label = label_list
        self.barcode_count = label_list.count('barcode')

    def output_objects(self, output_directory='path'):
        for sample in self.sample_list:
            object_list = []
            for count in range(self.barcode_count):
                object_list.append(open(output_directory + sample + '_' +
                                        str(count + 1) + '.fastq', 'w'))
            self.output_dict[sample] = object_list
        object_list = []
        for count in range(self.barcode_count):
            object_list.append(open(output_directory + 'unmatched' + '_' +
                                    str(count + 1) + '.fastq', 'w'))
        self.output_dict['unmatched'] = object_list

    def iterate_through_gseq(self):
        # transpose iterator list
        self.file_list = list(map(list, zip(*self.file_list)))
        for files in self.file_list:
            iterator = MultipleSequencingFileIterator(*files, directory=self.directory)
            barcode_indexs = duplicates(self.file_label, 'barcode')
            read_indexs = duplicates(self.file_label, 'read')
            for count, line in enumerate(iterator.iterator_zip()):
                self.reads += 1
                combined_filter = ''.join([qual[-1] for qual in line])
                if '0' not in combined_filter:
                    self.reads_pass_filter += 1
                    try:
                        key1 = self.barcode_1[line[barcode_indexs[0]][8]]
                    except KeyError:
                        key1 = None
                    try:
                        key2 = self.barcode_2[line[barcode_indexs[1]][8]]
                    except KeyError:
                        key2 = None
                    if key1 and key2:
                        sample_id = 'key' + str(key1) + ' key' + str(key2)
                        try:
                            sample = self.sample_key[sample_id]
                        except KeyError:
                            self.unmatched_read += 1
                            sample = 'unmatched'
                        out = self.output_dict[sample]
                        out[0].write(qseq_fastq_conversion(line[read_indexs[0]]))
                        out[1].write(qseq_fastq_conversion(line[read_indexs[1]]))
        for sample in self.output_dict.values():
            for out_object in sample:
                out_object.close()
