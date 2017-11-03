#!/usr/bin/env python3
from os import listdir

class DirectoryFileParser:
    """Opens Illumina qseq directory and processes qseq files, outputs samples fastq files"""

    def __init__(self, directory = 'path', file_description = []):
        self.directory = directory
        self.file_description = file_description

    def get_files(self):
        """Link to directory, pull  list of files, combine files with same unique id. Function only works with files
        that have an integer as a unique id
        -----------------------------------------------------
        self.directory: path to illumina sequencing lane directory
        self.file_description: list of lists containing sequencing file prefix and suffix, sequencing file designated
        by prefix*suffix
        returns; sorted list of relevant files names in directory"""
        files = []
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
            files.append([x for x, y in sorted(zip(sample_names[count], sorting_key[count]))])

        return files