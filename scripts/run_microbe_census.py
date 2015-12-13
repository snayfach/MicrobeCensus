#!/usr/bin/python

# MicrobeCensus - estimation of average genome size from shotgun sequence data
# Copyright (C) 2013-2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '1.0.5'

import argparse
from microbe_census import microbe_census

if __name__ == '__main__':

	# parse arguments
	parser = argparse.ArgumentParser(description='Estimate average genome size from metagenomic data.')

	io = parser.add_argument_group('Input/Output')
	parser.add_argument('seqfiles', type=str,
						help="""
							path to input metagenome(s);
							for paired-end metagenomes use commas to specify each file (ex: read_1.fq.gz,read_2.fq.gz);
							can be FASTQ/FASTA; can be gzip (.gz) or bzip (.bz2) compressed'
							""")
	io.add_argument('outfile', type=str, help="path to output file containing results")
	
	speed = parser.add_argument_group('Pipeline throughput (optional)')
	speed.add_argument('-n', dest='nreads', type=int, default=2000000,
						help="number of reads to sample from seqfile and use for AGS estimation. to use all reads set to 100000000. (default = 2000000)")
	speed.add_argument('-t', dest='threads', type=int, default=1,
						help="number of threads to use for database search (default = 1)")

	type = parser.add_argument_group('File type (optional)')
	type.add_argument('-f', dest='file_type', type=str,
						help="file type (default = autodetect)",
						choices=['fasta', 'fastq'])
	type.add_argument('-c', dest='fastq_format', type=str,
						help="quality score encoding (default = autodetect)",
						choices=['fastq-sanger', 'fastq-solexa', 'fastq-illumina'])

	qc = parser.add_argument_group('Quality control (optional)')
	qc.add_argument('-l', dest='read_length', type=int,
						help="all reads trimmed to this length; reads shorter than this discarded (default = median read length)",
						choices=[50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500])
	qc.add_argument('-q', dest='min_quality', type=int, default=-5,
						help="minimum base-level PHRED quality score (default = -5; no filtering)")
	qc.add_argument('-m', dest='mean_quality', type=int, default=-5,
						help="minimum read-level PHRED quality score (default = -5; no filtering)")
	qc.add_argument('-d', dest='filter_dups', action='store_true', default=False,
						help="filter duplicate reads (default = False)")
	qc.add_argument('-u', dest='max_unknown', type=int, default=100,
						help="max percent of unknown bases per read (default = 100 percent; no filtering)")

	parser.add_argument('-v', dest='verbose', action='store_true', default=False,
						help="print program\'s progress to stdout (default = False)")
	parser.add_argument('-V', '--version', action='version',
						version='MicrobeCensus (version %s)' % __version__)
	
	args = vars(parser.parse_args())
	args['seqfiles'] = args['seqfiles'].split(',')

	# run pipeline
	est_ags, args = microbe_census.run_pipeline(args)
	count_bases = microbe_census.count_bases(args['seqfiles'])

	# write output
	microbe_census.report_results(args, est_ags, count_bases)
