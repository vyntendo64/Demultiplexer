#!/usr/bin/env python3

import unittest
import Demultiplex
import subprocess

test_single_index_demultiplex = Demultiplex.Demuliplex(files = [
    {
        'path': 'tests/test_qseq/1_test.^.qseq.txt', 
        'action': 'read'
    }, 
    {
        'path': 'tests/test_qseq/2_test.^.qseq.txt', 
        'action': 'barcode'
    }],
    sample_key_path = 'tests/test_sample_files/single_index_test.txt',
    primary_barcodes_path = 'tests/test_sample_files/N700_nextera_barcodes.txt')

# print('test_single_index_demultiplex')
# print(vars(test_single_index_demultiplex))

test_single_index_demultiplex.run()
# test_single_index_demultiplex.process_barcodes()
# test_single_index_demultiplex.set_action()
# test_single_index_demultiplex.get_sample_labels()
# test_single_index_demultiplex.output_objects()
# test_single_index_demultiplex.iterate_through_qseq()

# test_single_index_demultiplex.run()


test_dual_index_demultiplex = Demultiplex.Demuliplex(files = [
    {
        'path': '1_test.^.qseq.txt',
        'action': 'read'
    },
    {
        'path': '2_test.^.qseq.txt',
        'action': 'barcode'
    },
    {
        'path': '3_test.^.qseq.txt',
        'action': 'barcode'
    },
    {
        'path': '4_test.^.qseq.txt',
        'action': 'read'
    }],
    sample_key_path = 'tests/test_sample_files/dual_index_test.txt',
    primary_barcodes_path = 'tests/test_sample_files/N700_nextera_barcodes.txt',
    secondary_barcodes_path = 'tests/test_sample_files/N500_nextera_barcodes.txt')

# test_dual_index_demultiplex.run()


# test_dual_index_demultiplex.run()


class TestDemultiplex(unittest.TestCase):

    def setUp(self):
        pass

    # # def test_reverse_complement(self):
    # #     self.assertEqual(DemultiplexClass.reverse_complement('GGCTATA'), 'TATAGCC')

    # def test_single_filter_pass(self):
    #     self.assertEqual(test_single_index_demultiplex.reads_pass_filter, 19614)

    # def test_single_unmatched_reads(self):
    #     self.assertEqual(test_single_index_demultiplex.unmatched_read, 2774)

    # def test_single_index_reads(self):
    #     self.assertEqual(test_single_index_demultiplex.indexed_reads, 16840)

    # def test_single_samples(self):
    #     x = test_single_index_demultiplex.file_list[0][0].replace('1_test', '')
    #     y = test_single_index_demultiplex.file_list[0][1].replace('2_test', '')
    #     self.assertEqual(x, y)

    # def test_dual_filter_pass(self):
    #     self.assertEqual(test_dual_index_demultiplex.reads_pass_filter, 19264)

    # def test_total_reads(self):
    #     self.assertEqual(test_single_index_demultiplex.reads, test_dual_index_demultiplex.reads)

    # def test_dual_unmatched_reads(self):
    #     self.assertEqual(test_dual_index_demultiplex.unmatched_read, 9420)

    # def test_dual_index_reads(self):
    #     self.assertEqual(test_dual_index_demultiplex.indexed_reads, 9844)

    # def test_dual_samples(self):
    #     x = test_dual_index_demultiplex.file_list[1][2].replace('3_test', '')
    #     y = test_dual_index_demultiplex.file_list[1][3].replace('4_test', '')
    #     self.assertEqual(x, y)

    # def test_command_line_single_index(self):
    #     parser_open = subprocess.run(
    #     [
    #         'python3', 
    #         'DemultiplexRun.py', 
    #         '-D', 
    #         'tests/test_qseq/', 
    #         '-S',
    #         'tests/test_sample_files/single_index_test.txt', 
    #         '-B1',
    #         'tests/test_sample_files/N700_nextera_barcodes.txt', 
    #         '-L', 
    #         'rb', 
    #         '-O', 
    #         'tests/test_output/',
    #         '-I', 
    #         '1_test.^.qseq.txt', 
    #         '2_test.^.qseq.txt'
    #     ],
    #     stdout = subprocess.PIPE)

    #     output = parser_open.stdout

    #     output_categories = (output.decode()).split('\n')
    #     filter = int(output_categories[2].split(':')[1])
    #     indexed = int(output_categories[3].split(':')[1])
    #     unmatched = int(output_categories[4].split(':')[1])

    #     self.assertEqual(filter, 19614)
    #     self.assertEqual(indexed, 16840)
    #     self.assertEqual(unmatched, 2774)

    # def test_command_line_dual_index(self):
    #     parser_open = subprocess.run(
    #     [
    #         'python3', 
    #         'DemultiplexRun.py', 
    #         '-D', 
    #         'tests/test_qseq/', 
    #         '-S',
    #         'tests/test_sample_files/dual_index_test.txt', 
    #         '-B1',
    #         'tests/test_sample_files/N700_nextera_barcodes.txt', 
    #         '-B2',
    #         'tests/test_sample_files/N500_nextera_barcodes.txt', 
    #         '-L', 
    #         'rbbr', 
    #         '-O',
    #         'tests/test_output/', 
    #         '-I', 
    #         '1_test.^.qseq.txt', 
    #         '2_test.^.qseq.txt',
    #         '3_test.^.qseq.txt', 
    #         '4_test.^.qseq.txt'
    #     ],
    #     stdout = subprocess.PIPE)

    #     output = parser_open.stdout
    #     print(output)
    #     output_categories = (output.decode()).split('\n')
    #     filter = int(output_categories[2].split(':')[1])
    #     indexed = int(output_categories[3].split(':')[1])
    #     unmatched = int(output_categories[4].split(':')[1])
    #     self.assertEqual(filter, 19264)
    #     self.assertEqual(indexed, 9420)
    #     self.assertEqual(unmatched, 9844)


if __name__ == '__main__':
    unittest.main()
