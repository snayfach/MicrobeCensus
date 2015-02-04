# MicrobeCensus - estimation of average genome size from shotgun sequence data
# Copyright (C) 2013-2014 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

#!/usr/bin/python

# import libraries
try:
    import sys
except Exception:
    print 'Module "sys" not installed'; exit()

try:
    import os
except Exception:
    print 'Module "os" not installed'; exit()
    
try:
    import optparse
except Exception:
    print 'Module "optparse" not installed'; exit()

try:
    import datetime
except Exception:
    print 'Module "datetime" not installed'; exit()

from sim_functions import *

#   parse options
parser = optparse.OptionParser(usage = "Usage: python seq_sim.py [-options] <seqfile> <outfile>")
parser.add_option("-l", "--read_length", dest="read_length", type="int", help="read length")
parser.add_option("-i", "--insert_dist", dest="insert_dist", type="int", help="insert distance")
parser.add_option("-n", "--n_reads", dest="n_reads", type="int", help="number of reads")
parser.add_option("-c", "--cov", dest="coverage", type="float", help="library coverage")
parser.add_option("-p", "--pe", dest="paired_end", default=False, action="store_true", help="paired end reads")
parser.add_option("-e", "--error_model", dest="error_model", default=None, help="type of sequencing error [illumina, uniform]")
parser.add_option("-r", "--error_rate", dest="error_rate", type="float", default=None, help="error rate (only valid for uniform error)")

(options, args) = parser.parse_args()
try:
	p_genome = args[0]
	f_out = open(args[1], 'w')
	read_length = options.read_length
	insert_dist = options.insert_dist
	n_reads = options.n_reads
	coverage = options.coverage
	paired_end = options.paired_end
	error_model = options.error_model
	error_rate = options.error_rate
except Exception, e:
    print "Incorrect number of command line arguments."
    print "\nUsage: python seq_sim.py [-options] <seqfile> <outfile>"
    print "For all options enter: python seq_sim.py -h"
    sys.exit()

# check for argument validity
if n_reads is None and coverage is None:
    sys.exit("Must specify library size using either --n_reads or --cov")
if read_length is None:
    sys.exit("Must specify read length using -l or --read_length")
if paired_end and not insert_dist:
    sys.exit("Must specify insert distance if simulating paired end library using -i or --insert_dist")

################################
# Main

# read in bias models
main_dir       = os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[0])))
src_dir        = os.path.join(main_dir, 'src')
data_dir       = os.path.join(main_dir, 'data')

# read in genome
genome = read_in_genome(p_genome)
scaffolds, weights, genome_size = parse_genome(genome)

# print library details
coverage = round((n_reads * read_length)/float(genome_size),4) if coverage is None else coverage
n_reads  = round(coverage * genome_size/float(read_length)) if n_reads is None else n_reads
insert_dist = read_length if not paired_end and not insert_dist else insert_dist
sys.stderr.write("DNA shotgun library:\n")  
sys.stderr.write("   Reference          = " + os.path.basename(p_genome) + "\n")
sys.stderr.write("   Read length        = " + str(read_length) + "\n")
sys.stderr.write("   Paired end reads   = " + str(paired_end) + "\n")
sys.stderr.write("   Insert distance    = " + (str(insert_dist) if paired_end else 'NA') + "\n")
sys.stderr.write("   Number of reads    = " + str(n_reads) + "\n")
sys.stderr.write("   Library coverage   = " + str(coverage) + "\n")
sys.stderr.write("   Error model        = " + str(error_model) + "\n")
sys.stderr.write("   Error rate         = " + str(error_rate) + "\n")

# read ids
read_id = 0
while read_id < n_reads:

	# pick scaffold/chromosome
	scaffold_seq = genome[weighted_choice(scaffolds, weights)]
	
	# pick fragment
	frag_start = random.randint(0, len(scaffold_seq) - 1)
	frag_stop    = frag_start + insert_dist
	frag_seq     = scaffold_seq[ frag_start : frag_stop ]

	# sequencing biases
	if len(frag_seq) < insert_dist: # fragment length
		continue

	# read is paired
	if paired_end:
		# sequence read
		pe1_seq = frag_seq[0:read_length]
		pe2_seq = revcomp(frag_seq[-read_length:])
		# generate errors
		if error_model:
			pe1_seq = mutate_read(pe1_seq, error_model, error_rate)
			pe2_seq = mutate_read(pe2_seq, error_model, error_rate)
		# write read to outfile
		f_out.write('>'+str(read_id)+'/1'+'\n'+str(pe1_seq)+'\n')
		f_out.write('>'+str(read_id)+'/2'+'\n'+str(pe2_seq)+'\n')

	# read is single-ended
	else:
		# sequence read
		pe1_seq = frag_seq[0:read_length]
		# generate errors
		if error_model:
			pe1_seq = mutate_read(pe1_seq, error_model, error_rate)
		# write read to outfile
		f_out.write('>'+str(read_id)+'\n'+str(pe1_seq)+'\n')

	# increment read id
	read_id += 1




