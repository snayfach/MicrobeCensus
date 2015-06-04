# MicrobeCensus - estimation of average genome size from shotgun sequence data
# Copyright (C) 2013-2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

"""
Usage:
	python sim_reads.py <threads> <read_length> <coverage>
	
	<threads>: number of threads to use
	<read_length>: read length to use for simulated shotgun libraries
	<coverage>: depth of coverage for each library

Description: simulate shotgun sequence libraries
	do for each genome listed in ../data/genomes
	genomes should be formatted as "genome_name".fna.gz
	simulated libraries are written to ../data/reads/read_length
	simulated libraries are formatted as "genome_name"-reads.fa
"""

# Libraries
import sys, os
from training import *

# Arguments
threads = int(sys.argv[1])
read_length = sys.argv[2]
coverage = sys.argv[3]

# Filepaths to input/intermediate/output data
package_dir  = os.path.dirname(os.path.realpath(sys.argv[0]))
input_dir    = os.path.join(package_dir, 'input')
inter_dir    = os.path.join(package_dir, 'intermediate')
out_dir      = os.path.join(package_dir, 'output')
genomes_dir  = os.path.join(input_dir, 'genomes')
reads_dir    = os.path.join(inter_dir, 'reads')
lib_dir      = os.path.join(reads_dir, read_length)
readlen_path = os.path.join(out_dir, 'read_len.map')

# Check if input directories exist
for dir in [input_dir, inter_dir, out_dir, genomes_dir, reads_dir, lib_dir]:
	if not os.path.isdir(dir):
		os.mkdir(dir)

# Build commands to run in parallel
command   = "python %(package_dir)s/seq_sim.py --read_length=%(read_length)s --cov=%(coverage)s %(genome_path)s %(lib_path)s"
args_list = []
for genome_file in os.listdir(genomes_dir):
	if genome_file[-7:] != '.fna.gz': # check genome file extension
		sys.exit('Genome files must be formatted: "genome_name".fna.gz')
	else:
		genome_name = genome_file[0:-7]
		genome_path = os.path.join(genomes_dir, genome_file)
		lib_path = os.path.join(lib_dir, genome_name + '-reads.fa')
		args = {'package_dir':package_dir, 'genome_path':genome_path, 'lib_path':lib_path, 'coverage':coverage, 'read_length':read_length}
		args_list.append(args)

# Print arguments
print "Threads to use: %s" % threads
print "Read length: %s" % read_length
print "Genome coverage: %s" % coverage
print "Number of libraries: %s" % len(args_list)
print "Output directory: %s" % reads_dir

# Simulate libraries in parallel
parallel_subprocess(command, args_list, threads)

# Record read length
f_out = open(readlen_path, 'a')
f_out.write(str(read_length) + '\n')







