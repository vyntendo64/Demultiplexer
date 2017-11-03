#!/usr/bin/env python3

class Barcode:
    def __init__(self, barcode, number):
        self.barcode = barcode
        self.number = number
        self.mismatch_list = ['A', 'T', 'G', 'C', '.']
        self.complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


    def reverse(self):
        """Simple reverse complement function used to initialize a barcode dictionary, (original sequence and reverse
        complement both link to same hash in dictionary).
        -----------------------------------------------------
        string='string': string must be composed of BP ATGC
        returns; reverse complement string"""
        # reverse string
        reversed = []
        self.barcode = self.barcode[::-1]
        # complementary bp lookup dictionary
        bases = list(self.barcode)
        # iterate over string list
        for i in bases:
            reversed.append(self.complement[i])

        # return joined string
        self.barcode = ''.join(reversed)

        return self.barcode

    def get_possible_mismatches(self, mismatch):
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

        # print('mismatch list')
        # print(self.mismatch_list)

        all_barcodes = [self.barcode]

        for possible_mismatch in self.mismatch_list:
            # print('possible mismatch')
            # print(possible_mismatch)
            # loop over every position in barcode
            for index in range(len(self.barcode)):
                # barcode to list
                barcode_list = list(self.barcode)
                # set character in to mismatch
                barcode_list[index] = possible_mismatch
                # if mismatched barcode not already in list then add the barcode to the list
                if ''.join(barcode_list) not in all_barcodes:
                    all_barcodes.append(''.join(barcode_list))


        if mismatch > 1:
            # initialize additional loops based on mismatch number
            for count in range(mismatch - 1):
                # repeat loop
                for mismatched_barcode in list(all_barcodes):
                    for possible_mismatch in self.mismatch_list:
                        for character in range(len(mismatched_barcode)):
                            barcode_list = list(mismatched_barcode)
                            barcode_list[character] = possible_mismatch
                            if ''.join(barcode_list) not in all_barcodes:
                                all_barcodes.append(''.join(barcode_list))

        return list(set(all_barcodes))