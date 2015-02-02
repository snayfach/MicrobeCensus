#!/usr/bin/python

import argparse
from microbe_census import microbe_census

if __name__ == '__main__':

	# parse arguments
	parser = argparse.ArgumentParser(description='Estimate average genome size from metagenomic data.')
	parser.add_argument('seqfile', type=str,
						help='path to input metagenome')
	parser.add_argument('outfile', type=str,
						help='path to output file')
	parser.add_argument('-n', dest='nreads', type=int, default=1000000,
						help='number of reads to use for AGS estimation (default = 1e6)')
	parser.add_argument('-l', dest='read_length', type=int,
						help='trim reads to this length (default = median read length)',
						choices=[50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500])
	parser.add_argument('-f', dest='file_type', type=str,
						help='file type (default = autodetect)',
						choices=['fasta', 'fastq'])
	parser.add_argument('-c', dest='fastq_format', type=str,
						help='quality score encoding (default = autodetect)',
						choices=['fastq-sanger', 'fastq-solexa', 'fastq-illumina'])
	parser.add_argument('-t', dest='threads', type=int, default=1,
						help='number of threads to use for database search (default = 1)')
	parser.add_argument('-q', dest='min_quality', type=int, default=-5,
						help='minimum base-level PHRED quality score (default = -5; no filtering)')
	parser.add_argument('-m', dest='mean_quality', type=int, default=-5,
						help='minimum read-level PHRED quality score (default = -5; no filtering)')
	parser.add_argument('-d', dest='filter_dups', action='store_true', default=False,
						help='filter duplicate reads (default = False)')
	parser.add_argument('-u', dest='max_unknown', type=int, default=100,
						help='max percent of unknown bases (default = 100 percent; no filtering)')
	parser.add_argument('-k', dest='keep_tmp', action='store_true', default=False,
						help='keep temporary files (default = False)')
	parser.add_argument('-s', dest='quiet', action='store_true', default=False,
						help='suppress printing program\'s progress to stdout (default = False)')
	args = vars(parser.parse_args())

	# run pipeline
	est_ags, args = microbe_census.run_pipeline(args)

	# write output
	microbe_census.report_results(args, est_ags)