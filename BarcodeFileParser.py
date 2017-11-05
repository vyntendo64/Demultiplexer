#!/usr/bin/env python3
from Barcode import Barcode

class BarcodeFileParser:
    def __init__(self, path = None, 
        mismatch = 1,
        mismatch_list = ['A', 'T', 'G', 'C', '.']):

        self.path = path
        self.mismatch = 0 if mismatch <= 0 else mismatch
        self.mismatch_list = mismatch_list

    def set_mismatch_list(self, mismatch_list = []):
        self.mismatch_list = mismatch_list

    def get_barcodes(self):
        barcodes = {}

        if (self.path == None):
            return barcodes

        for index, line in enumerate(open(self.path)):

            barcode = self.read_barcode_from_line(line = line)
            
            mismatches = self.get_possible_mismatches(barcode = barcode.get())
            barcodes.update(self.get_barcodes_from_mismatches(mismatches = mismatches, 
                barcode_number = barcode.get_number()))

            mismatches = self.get_possible_mismatches(barcode = barcode.reverse())
            barcodes.update(self.get_barcodes_from_mismatches(mismatches = mismatches, 
                barcode_number = barcode.get_number()))

        return barcodes

    def read_barcode_from_line(self, line = None):
        parsed = (line.replace('\n', '')).split('\t')

        return Barcode(parsed[0], int(parsed[1]))

    def get_possible_mismatches(self, barcode):
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
        # calculate first set of mismatches
        # loop over possible mismatches
        
        mismatched_barcodes = [barcode]
        # barcode to list

        for possible_mismatch in self.mismatch_list:
            # print('possible mismatch')
            # print(possible_mismatch)
            # loop over every position in barcode
            for index in range(len(barcode)):        
                barcode_list = list(barcode)
                # set character in to mismatch
                barcode_list[index] = possible_mismatch
                # if mismatched barcode not already in list then add the barcode to the list
                if ''.join(barcode_list) not in mismatched_barcodes:
                    mismatched_barcodes.append(''.join(barcode_list))

        if self.mismatch > 1:
            # initialize additional loops based on mismatch number
            for count in range(self.mismatch - 1):
                # repeat loop
                for mismatched_barcode in list(mismatched_barcodes):
                    for possible_mismatch in self.mismatch_list:
                        for character in range(len(mismatched_barcode)):

                            barcode_list = list(mismatched_barcode)
                            barcode_list[character] = possible_mismatch
                            if ''.join(barcode_list) not in mismatched_barcodes:
                                mismatched_barcodes.append(''.join(barcode_list))

        return list(set(mismatched_barcodes))

    def get_barcodes_from_mismatches(self, mismatches, barcode_number):
        barcodes = {}
        collision_sequences = []
        for possible_barcode in mismatches:
            if possible_barcode not in collision_sequences:
                if possible_barcode not in list(barcodes.keys()):
                    barcodes[possible_barcode] = barcode_number
                else:
                    collision_sequences.append(possible_barcode)
                    del barcodes[possible_barcode]
        return barcodes
