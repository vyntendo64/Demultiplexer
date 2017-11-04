#!/usr/bin/env python3

import gzip
import io

class MultipleSequencingFileIterator:
    def __init__(self, 
        files = [], 
        directory = 'path', 
        gnu_zipped = False):
        """Initiate iteration object, yield line in gseq files
         -----------------------------------------------------
         *args='path_to_gesq': returns an iterator object for paired sequencing files
         gnu_zipped=False: if gnu_zipped=True will process files with python gzip library, increases processing time
         unix gunzip is faster"""

        self.iter_list = []

        for file in files:
            print(file)
            print('file')
            self.iter_list.append(self.iteration_call(file["path"]))

    def iteration_call(iter_file):
        if gnu_zipped:
            # wrap encoded line in buffer to stream text file
            with io.TextIOWrapper(io.BufferedReader(gzip.open(iter_file, 'rb'))) as seq:
                for line in seq:
                    # yield a line split by tabs and stripped of line identifier, '\n'
                    yield ((line.replace('\n', '')).split('\t'))
        else:
            with open(iter_file) as seq:
                for line in seq:
                    # yield a line split by tabs and stripped of line identifier, '\n'
                    yield ((line.replace('\n', '')).split('\t'))

    # zip files together to iterate over all of the files at once and yield one line at a time for looping
    def iterator_zip(self):
        for line in zip(*self.iter_list):
            yield line
