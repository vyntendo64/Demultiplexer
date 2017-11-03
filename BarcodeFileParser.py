#!/usr/bin/env python3
from Barcode import Barcode

class BarcodeFileParser:
    """Open multiple fastq or gseq files together and iterate over them as a group"""

    def __init__(self, path = None, mismatch = 1):
        self.path = path
        self.mismatch = 0 if mismatch <= 0 else mismatch

        self.set_up()

    def set_up(self):
        self.collision_sequences = []
        self.barcodes = {}

    def get_barcodes(self):
        # print(self.path)
        if (self.path == None):
            return self.barcodes

        for index, line in enumerate(open(self.path)):

            barcode = self.read_barcode_from_line(line)
            self.set_barcodes_from_list(barcode.get_possible_mismatches(self.mismatch), barcode)

            # repeat for barcode reverse complement
            barcode.reverse()
            self.set_barcodes_from_list(barcode.get_possible_mismatches(self.mismatch), barcode)

        return self.barcodes

    def read_barcode_from_line(self, line):
        parsed = (line.replace('\n', '')).split('\t')

        return Barcode(parsed[0], int(parsed[1]))

    def set_barcodes_from_list(self, listing, barcode):
        for possible_barcode in listing:
            if possible_barcode not in self.collision_sequences:
                if possible_barcode not in list(self.barcodes.keys()):
                    self.barcodes[possible_barcode] = barcode.number
                else:
                    self.collision_sequences.append(possible_barcode)
                    del self.barcodes[possible_barcode]