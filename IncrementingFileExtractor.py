#!/usr/bin/env python3
import os 
import sys

class IncrementingFileExtractor:
    def __init__(self, files = []):
        self.files = files
        self.break_characters = ['.^.']

    def get_files(self):
        files = []

        for file in self.files:
            base = os.path.basename(file["path"])

            for break_character in self.break_characters:

                if break_character in base:
                    incremented_files = self.get_all_incrementing_files(file)                
                    files.append({'path': incremented_files, 'action': file["action"]})

                else:
                    files.append({'path': file["path"], 'action': file['action']})
        return files

    def get_all_incrementing_files(self, file):
        files = []
        path = file["path"]
        base = os.path.basename(path)
        directory = os.path.dirname(path)
        file_description = base.split('^')
        incrementor = 1
        is_file = True

        while(is_file):
            possible_file = file_description[0] + str(incrementor) + file_description[1]

            if(os.path.isfile(directory + '/' + possible_file)):
                files.append(directory + '/' + possible_file)
                incrementor = incrementor + 1
                
            else:
                is_file = False

        return files

    def get_read_count(self):
        read_count = 0;

        for file in self.files:
            if(file["action"] == 'read'):
                read_count = read_count + 1

        return read_count

    def get_barcode_count(self):
        barcode_count = 0;

        for file in self.files:
            if(file["action"] == 'barcode'):
                barcode_count = barcode_count + 1

        return barcode_count