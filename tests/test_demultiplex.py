#!/usr/bin/env python3

import unittest
import DemultiplexClass

test_single_index_demultiplex = DemultiplexClass.Demuliplex('1_test.*.qseq.txt', '2_test.*.qseq.txt',
                                                            directory='test_qseq/',
                                                            sample_key='test_sample_files/single_index_test.txt',
                                                            barcode_1='test_sample_files/N700_nextera_barcodes.txt',
                                                            file_label='rb',
                                                            )
test_single_index_demultiplex.get_directory_lists()
test_single_index_demultiplex.process_barcodes()
test_single_index_demultiplex.process_file_label()
test_single_index_demultiplex.get_sample_labels()
test_single_index_demultiplex.output_objects(output_directory='test_output/')
test_single_index_demultiplex.iterate_through_qseq()


test_dual_index_demultiplex = DemultiplexClass.Demuliplex('1_test.*.qseq.txt', '2_test.*.qseq.txt',
                                                          '3_test.*.qseq.txt', '4_test.*.qseq.txt',
                                                          directory='test_qseq/',
                                                          sample_key='test_sample_files/dual_index_test.txt',
                                                          barcode_1='test_sample_files/N700_nextera_barcodes.txt',
                                                          barcode_2='test_sample_files/N500_nextera_barcodes.txt',
                                                          file_label='rbbr'
                                                          )
test_dual_index_demultiplex.get_directory_lists()
test_dual_index_demultiplex.process_barcodes()
test_dual_index_demultiplex.process_file_label()
test_dual_index_demultiplex.get_sample_labels()
test_dual_index_demultiplex.output_objects(output_directory='test_output/')
test_dual_index_demultiplex.iterate_through_qseq()


class TestDemultiplex(unittest.TestCase):

    def setUp(self):
        pass

    def test_reverse_complement(self):
        self.assertEqual(DemultiplexClass.reverse_complement('GGCTATA'), 'TATAGCC')

    def test_single_filter_pass(self):
        self.assertEqual(test_single_index_demultiplex.reads_pass_filter, 19614)

    def test_single_unmatched_reads(self):
        self.assertEqual(test_single_index_demultiplex.unmatched_read, 2774)

    def test_single_index_reads(self):
        self.assertEqual(test_single_index_demultiplex.indexed_reads, 16840)

    def test_single_samples(self):
        x = test_single_index_demultiplex.file_list[0][0].replace('1_test', '')
        y = test_single_index_demultiplex.file_list[0][1].replace('2_test', '')
        self.assertEqual(x, y)

    def test_dual_filter_pass(self):
        self.assertEqual(test_dual_index_demultiplex.reads_pass_filter, 19264)

    def test_total_reads(self):
        self.assertEqual(test_single_index_demultiplex.reads, test_dual_index_demultiplex.reads)

    def test_dual_unmatched_reads(self):
        self.assertEqual(test_dual_index_demultiplex.unmatched_read, 9420)

    def test_dual_index_reads(self):
        self.assertEqual(test_dual_index_demultiplex.indexed_reads, 9844)

    def test_dual_samples(self):
        x = test_dual_index_demultiplex.file_list[1][2].replace('3_test', '')
        y = test_dual_index_demultiplex.file_list[1][3].replace('4_test', '')
        self.assertEqual(x, y)

if __name__ == '__main__':
    unittest.main()
