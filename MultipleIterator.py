#!/usr/bin/env python3

import gzip
import io

class MultipleSequencingFileIterator:
    def __init__(self, 
        files = [],
        gnu_zipped = False):

        self.files = files
        self.gnu_zipped = gnu_zipped
        self.build = []
        self.ordered_paths = self.get_ordered_paths(files)
        self.sequence()

    def get_ordered_paths(self, files):
        ordered_paths = []

        for file in files:
            ordered_paths.append(file['path'])

        return list(map(list, zip(*ordered_paths)))

    def sequence(self):
        for path_no, ordered_path in enumerate(self.ordered_paths):
            rows = []

            for path in ordered_path:
                rows.append(self.get_row_from_file_path(path))

            self.build.append(rows)

    def get_row_from_file_path(self, path):
        if self.gnu_zipped:
            # wrap encoded row in buffer to stream text file
            with io.TextIOWrapper(io.BufferedReader(gzip.open(path, 'rb'))) as file:
                for row in file:
                    # yield a row split by tabs and stripped of row identifier, '\n'
                    yield ((row.replace('\n', '')).split('\t'))
        else:
            with open(path) as file:
                for row in file:
                    # yield a row split by tabs and stripped of row identifier, '\n'
                    yield ((row.replace('\n', '')).split('\t'))

    def get_next_row(self):
        for file in self.build:
            for row in zip(*file):
                yield row

    def get_barcode_indexes(self):
        barcode_indexes = []

        for index, file in enumerate(self.files):

            if file['action'] == 'barcode':
                barcode_indexes.append(index)

        return barcode_indexes

    def get_read_indexes(self):
        read_indexes = []

        for index, file in enumerate(self.files):

            if file['action'] == 'read':
                read_indexes.append(index)

        return read_indexes
