# MicrobeCensus - estimation of average genome size from shotgun sequence data
# Copyright (C) 2013-2014 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

"""
Usage:
	python class_reads.py <threads>
	
	<threads>: number of threads to use

Description: classify reads into marker gene families across range of mapping parameters
	do for all m8 files listed in ../data/search
	m8 files should be formatted as "genome_name".m8
	classification results are written to ../data/hits/read_length
	classification results are formatted as "genome_name".hits

"""

# Libraries
import sys, os, Bio.SeqIO
from training import *

# Arguments
threads = int(sys.argv[1])

# Filepaths to input/intermediate/output data
main_dir      = os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[0])))
gene_fam_dir  = os.path.join(main_dir, 'training/input/gene_fams')
search_dir    = os.path.join(main_dir, 'training/intermediate/search')
hits_dir      = os.path.join(main_dir, 'training/intermediate/hits')
gene_fam_path = os.path.join(main_dir, 'training/output/gene_fam.map')
gene_len_path = os.path.join(main_dir, 'training/output/gene_len.map')

# create directories
if not os.path.isdir(hits_dir): os.mkdir(hits_dir)

# Read in lookups
gene2fam = {}
gene2len = {}
fams = set([])
for file in os.listdir(gene_fam_dir):
	if file[-7:] != '.faa.gz': sys.exit('Gene family multi-fasta files must have .faa.gz extension')
	else:
		fam = file[0:-7]
		fams.add(fam)
		for rec in Bio.SeqIO.parse(gzip.open(os.path.join(gene_fam_dir, file)), 'fasta'):
			gene2fam[str(rec.id)] = fam
			gene2len[str(rec.id)] = len(str(rec.seq))

# Setup cutoffs to test
aln_covs = [0.00, 0.25, 0.50, 0.75]
max_pids = [50, 60, 70, 80, 90, 100]
min_scores = drange(23, 50, 1)

# Build commands to run in parallel
args_list = []
for read_length in os.listdir(search_dir):
	in_dir  = os.path.join(search_dir, read_length)
	out_dir = os.path.join(hits_dir, read_length)
	if not os.path.isdir(out_dir): os.mkdir(out_dir) # create outdir if not exists
	for search_file in os.listdir(in_dir):
		if search_file[-3:] == '.m8': # check file extension
			genome_name = search_file[0:-3]
			search_path = os.path.join(in_dir, search_file)
			hits_path = os.path.join(out_dir, genome_name + '.hits')
			args = {'p_in':search_path, 'p_out':hits_path, 'aln_covs':aln_covs, 'max_pids':max_pids, 'min_scores':min_scores, 'gene2len':gene2len, 'gene2fam':gene2fam, 'fams':fams, 'read_length':read_length}
			args_list.append(args)

# Print arguments
print "Threads to use: %s" % threads
print "Number of search results to classify: %s" % len(args_list)
print "Output directory: %s" % hits_dir

# Classify reads in parallel
parallel_function(classify_reads, args_list, threads)

# Write gene_fam.map and gene_len.map
f_out = open(gene_fam_path, 'w')
for gene, fam in gene2fam.iteritems():
	record = [str(gene), str(fam)]
	f_out.write('\t'.join(record)+'\n')
f_out = open(gene_len_path, 'w')
for gene, length in gene2len.iteritems():
	record = [str(gene), str(length)]
	f_out.write('\t'.join(record)+'\n')



