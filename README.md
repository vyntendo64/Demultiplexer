# Demultiplexer

Script demultiplexes Illumnina qseq lane files to sample .fastq files

## Usage

python3 DemultiplexRun -D directory -S sample_key -B1 barcode_1 -B2 barcode_2 -L file_labels -mismatch_number -O output_directory -I input_file_1 input_file_2 ...

## Inputs

- -D, /path/ to qseq directory
- -S, /path/sample_file.txt; file should be formatted as 'barcode tab sample_name' for single index and 'barcode tab barcode tab sample_name for dual indexes
- -B1, /path/barcode_1_file, line separated list of barcodes
- -B2, /path/barcode_2_file, line separated list of barcodes
- -L, string of r and b character to designate input files as barcode or read files, should be the same order as input files
- -M, number of barcode mismatches to consider, default = 2
- -O, path to output directory
- -Z, designate is inpute qseq files are gzipped, slows processing
- -I, qseq file prefix and suffix separated by *, ie. -I s_1_.*.qseq.txt s_2_.*.qseq.txt

## Examples

### Single Index Demultiplex

python3 DemultiplexRun -D ~/Demultiplexer/tests/test_qseq -S ~/Demultiplexer/tests/test_sample_files/single_index_test.txt -B1 ~/Demultiplexer/tests/test_sample_files/N700_nextera_bacrodes.txt -L 'rb' -M 1 -O ~/Demultiplexer/tests/test_output/ -I 1_test.*.qseq.txt 2_test.*.qseq.txt

### Dual Indes Demultiplex

python3 DemultiplexRun -D ~/Demultiplexer/tests/test_qseq -S ~/Demultiplexer/tests/test_sample_files/single_index_test.txt -B1 ~/Demultiplexer/tests/test_sample_files/N700_nextera_bacrodes.txt -B2 ~/Demultiplexer/tests/test_sample_files/N500_nextera_bacrodes.txt -L 'rbbr' -M 1 -O ~/Demultiplexer/tests/test_output/ -I 1_test.*.qseq.txt 2_test.*.qseq.txt 3_test.*.qseq.txt 4_test.*.qseq.txt
