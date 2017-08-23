#! /usr/env python3


class MultipleSequencingFileIterator:
    """Open multiple fastq or gseq files together and iterate over them as a group"""

    def __init__(self, *args):
        """Initiate iteration object, yield line in gseq files
         -----------------------------------------------------
         *args='path_to_gesq': returns an iterator object for paired sequencing files"""
        file_list = []
        for file in args:
            file_list.append(file)

        def iteration_call(iter_file):
            with open(iter_file) as seq:
                for line in seq:
                    # yield a line split by tabs and stripped of line identifier, '\n'
                    yield ((line.strip('\n')).split('\t'))

        self.iter_list = []
        # append iterator object to list
        for file in file_list:
            self.iter_list.append(iteration_call(file))

    def iterator_zip(self):
        for line in zip(*self.iter_list):
            yield line