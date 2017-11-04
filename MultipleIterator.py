#!/usr/bin/env python3

import gzip
import io

class MultipleSequencingFileIterator:
    def __init__(self, 
        files = []):
        self.gnu_zipped = False
        self.files = files
        self.ordered_paths = self.get_ordered_paths(files)
        self.build = []

    def populate(self):
        for ordered_path in self.ordered_paths:
            for path in ordered_path:
                self.build.append(self.get_row_from_file_path(path))

    def get_ordered_paths(self, files):
        ordered_paths = []
        for file in files:
            ordered_paths.append(file['path'])

        return list(map(list, zip(*ordered_paths)))
            

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
        for row in zip(*self.build):
            yield row
