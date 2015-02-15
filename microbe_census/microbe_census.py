# MicrobeCensus - estimation of average genome size from shotgun sequence data
# Copyright (C) 2013-2014 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

__version__ = '1.0.3'

#######################################################################################
#   IMPORT LIBRARIES

try: import sys
except: sys.exit("Could not import module 'sys'")
   
try: import os
except: sys.exit("Could not import module 'os'")
          
try: import argparse
except: sys.exit("Could not import module 'argparse'")

try: import platform
except: sys.exit("Could not import module 'platform'")

try: import gzip
except: sys.exit("Could not import module 'gzip'")

try: import subprocess
except: sys.exit("Could not import module 'subprocess'")

try: import bz2
except: sys.exit("Could not import module 'bz2'")

try: from Bio.SeqIO import parse
except: sys.exit("Could not import module 'Bio.SeqIO.parse'")

try: from numpy import median, mean, std
except: sys.exit("Could not import module 'numpy'")

try: from tempfile import mkstemp
except: sys.exit("Could not import module 'tempfile.mkstemp'")

try: import io
except: sys.exit("Could not import module 'io'")

try: import pkg_resources
except: sys.exit("Could not import module 'pkg_resources'")

#######################################################################################
#   FUNCTIONS

def mad(x, const = 1.48):
    """ Calculate the median absolute deviation (MAD) of x """
    return const * median([abs(i - median(x)) for i in x])

def open_file(inpath):
	""" Open input file for reading regardless of compression [gzip, bzip] or python version """
	ext = inpath.split('.')[-1]
	# Python2
	if sys.version_info[0] == 2:
		if ext == 'gz': return gzip.open(inpath)
		elif ext == 'bz2': return bz2.BZ2File(inpath)
		else: return open(inpath)
	# Python3
	elif sys.version_info[0] == 3:
		if ext == 'gz': return io.TextIOWrapper(gzip.open(inpath))
		elif ext == 'bz2': return bz2.BZ2File(inpath)
		else: return open(inpath)

def find_opt_pars(path_optpars, read_length):
    """ Read in optimal parameters for each family at given read length
        Returns optpars[fam_id] = [min_cov, max_aaid, min_score, aln_stat]
    """
    optpars = {}
    with open(path_optpars) as f_in:
        next(f_in)
        for line in f_in:
            fam_id, read_len, min_cov, max_aaid, min_score, aln_stat = line.rstrip().split()
            if read_len != str(read_length): continue
            else: optpars[fam_id] = {'min_cov':float(min_cov), 'max_aaid':float(max_aaid), 'min_score':float(min_score), 'aln_stat':aln_stat}
    return optpars

def read_dic(file, header, dtype, empty=False):
    """ Read in simple key value dictionary from file 
        Returns dic[field1] = field2
    """
    dic = {}
    f_in = gzip.open(file) if file.split('.')[-1] == 'gz' and file.split('.') > 1 else open(file)
    if header is True: next(f_in)
    for line in f_in:
        if empty is True:
            key = line.rstrip()
            value = None
        else:
            key, value = line.rstrip().split()
        dic[key] = float(value) if dtype == 'float' else int(value) if dtype == 'int' else value
    return dic

def read_list(file, header, dtype):
    """ Read in list from file """
    my_list = []
    f_in = gzip.open(file) if file.split('.')[-1] == 'gz' and file.split('.') > 1 else open(file)
    if header is True: next(f_in)
    for line in f_in:
        value = line.rstrip()
        my_list.append(float(value)) if dtype == 'float' else my_list.append(int(value)) if dtype == 'int' else my_list.append(value)
    return my_list

def check_os():
	""" Return os name """
	if platform.system() not in ['Linux', 'Darwin']:
		sys.exit("Operating system '%s' not supported" % platform.system())

def get_relative_paths():
	""" Fetch relative paths to data files """
	paths = {}
	paths['rapsearch'] = pkg_resources.resource_filename('microbe_census', 'bin/%s' % '_'.join(['rapsearch',platform.system(),'2.15']))
	paths['db'] = pkg_resources.resource_filename('microbe_census', 'data/rapdb_2.15')
	paths['fams'] = pkg_resources.resource_filename('microbe_census', 'data/gene_fam.map')
	paths['genelen'] = pkg_resources.resource_filename('microbe_census', 'data/gene_len.map')
	paths['params'] = pkg_resources.resource_filename('microbe_census', 'data/pars.map')
	paths['coeffs'] = pkg_resources.resource_filename('microbe_census', 'data/coefficients.map')
	paths['weights'] = pkg_resources.resource_filename('microbe_census', 'data/weights.map')
	paths['readlen'] = pkg_resources.resource_filename('microbe_census', 'data/read_len.map')
	paths['tempfile'] = mkstemp()[1]
	return paths
	
def check_paths(paths):
	""" Check that all relative paths exist """
	for my_path in paths.values():
		if os.path.isfile(my_path):
			continue
		elif os.path.isdir(my_path):
			continue
		else:
			sys.exit("Path to file/dir not found: %s" % my_path)

def auto_detect_read_length(seqfile, file_type):
	""" Find median read length from first 10K reads in seqfile """
	valid_lengths = [50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500]
	read_lengths = []
	for index, record in enumerate(parse(open_file(seqfile), file_type)):
		if index == 10000: break
		read_lengths.append(len(record.seq))
	median_read_length = int(median(read_lengths))
	if median_read_length < valid_lengths[0]:
		sys.exit("Median read length is %s. Cannot compute AGS using reads shorter than 50 bp." % median_read_length)
	for index, read_length in enumerate(valid_lengths):
		if read_length > median_read_length:
			return valid_lengths[index-1]
	return valid_lengths[-1]

def auto_detect_file_type(seqfile):
	""" Detect file type [fasta or fastq] of <p_reads> """
	for line in open_file(seqfile):
		if line[0] == '>': return 'fasta'
		elif line[0] == '@': return 'fastq'
		else: sys.exit("Filetype [fasta, fastq] of %s could not be recognized" % seqfile)

def auto_detect_fastq_format(seqfile):
	""" Use first 50,000 reads to detect quality encoding """
	max_reads = 50000
	formats = ['fastq-illumina', 'fastq-solexa', 'fastq-sanger']
	for format in formats:
		try:
			for index, rec in enumerate(parse(open_file(seqfile), format)):
				if index == max_reads: break
			return format
		except:
			pass
	sys.exit("Could not determine FASTQ quality encoding")

def impute_missing_args(args):
	""" Fill in missing arguments with defaults """
	if 'verbose' not in args:
		args['verbose'] = False
	
	if 'outfile' not in args:
		args['outfile'] = None
		
	if 'nreads' not in args:
		args['nreads'] = 1000000

	if 'threads' not in args:
		args['threads'] = 1

	if 'filter_dups' not in args:
		args['filter_dups'] = False
		
	if 'keep_tmp' not in args:
		args['keep_tmp'] = False
		
	if 'mean_quality' not in args:
		args['mean_quality'] = -5
		
	if 'min_quality' not in args:
		args['min_quality'] = -5

	if 'max_unknown' not in args:
		args['max_unknown'] = 100

	if 'file_type' not in args or args['file_type'] is None:
		args['file_type'] = auto_detect_file_type(args['seqfile'])

	if 'read_length' not in args or args['read_length'] is None:
		args['read_length'] = auto_detect_read_length(args['seqfile'], args['file_type'])

	if 'fastq_format' not in args:
		args['fastq_format'] = None

	if args['file_type'] == 'fastq' and args['fastq_format'] is None:
		args['fastq_format'] = auto_detect_fastq_format(args['seqfile'])
		
	if args['file_type'] == 'fastq':
		args['quality_type'] = 'solexa_quality' if args['fastq_format'] == 'fastq-solexa' else 'phred_quality'

def check_arguments(args):
	# QC options are not available for FASTA files
	if args['file_type'] == 'fasta' and any([args['min_quality'] > -5, args['mean_quality'] > -5, args['fastq_format'] is not None]):
		sys.exit("Quality filtering options are only available for FASTQ files")
	# check number of threads
	if args['threads'] < 1:
		sys.exit("Invalid number of threads: %s\nMust be a positive integer." % args['threads'])
	# check number of reads
	if args['nreads'] is not None and args['nreads'] < 1:
		sys.exit("Invalid number of reads: %s\nMust be a positive integer." % args['nreads'])

def print_copyright():
	# print out copyright information
	print ("\nMicrobeCensus - estimation of average genome size from shotgun sequence data")
	print ("version %s; github.com/snayfach/MicrobeCensus" % __version__)
	print ("Copyright (C) 2013-2014 Stephen Nayfach")
	print ("Freely distributed under the GNU General Public License (GPLv3)\n")

def print_parameters(args):
	# print out parameters
	print ("=============Parameters==============")
	print ("Input metagenome: %s" % args['seqfile'])
	print ("Output file: %s" % args['outfile'])
	print ("Reads trimmed to: %s bp" % args['read_length'])
	print ("Maximum reads sampled: %s" % args['nreads'])
	print ("Threads to use for db search: %s" % args['threads'])
	print ("File format: %s" % args['file_type'])
	print ("Quality score encoding: %s" % (args['fastq_format'] if args['file_type'] == 'fastq' else 'NA'))
	print ("Minimum base-level quality score: %s" % (args['min_quality'] if args['file_type'] == 'fastq' else 'NA'))
	print ("Minimum read-level quality score: %s" % (args['mean_quality'] if args['file_type'] == 'fastq' else 'NA'))
	print ("Maximum percent unknown bases/read: %s" % args['max_unknown'])
	print ("Filter duplicate reads: %s" % args['filter_dups'])
	print ("Keep temporary files: %s\n" % args['keep_tmp'])

def quality_filter(rec, args):
	""" Return true if read fails QC """
	# parse record
	sequence = rec.seq[0:args['read_length']]
	quality = rec.letter_annotations[args['quality_type']][0:args['read_length']] if args['file_type'] == 'fastq' else None
	# check percent unknown
	if 100 * sum([1 if b == 'N' else 0 for b in sequence]) / float(len(sequence)) > args['max_unknown']:
		return True
	# check mean quality
	elif quality and mean(quality) < args['mean_quality']:
		return True
	# check min quality
	elif quality and min(quality) < args['min_quality']:
		return True
	# read passed QC
	else:
		return False

def process_seqfile(args, paths):
	""" Sample high quality reads from seqfile """
	if args['verbose']:
		print ("====Estimating Average Genome Size====")
		print ("Sampling & trimming reads...")
	outfile = open(paths['tempfile'], 'w')
	# loop over sequences
	read_id, dups, too_short, low_qual = 0, 0, 0, 0
	seqs = set([])
	for rec in parse(open_file(args['seqfile']), args['fastq_format'] if args['file_type'] == 'fastq' else 'fasta'):
		# record sequence if enough high quality bases remain
		if len(rec.seq) < args['read_length']:
			too_short += 1; continue
		# check if sequence is a duplicate
		elif args['filter_dups'] and (str(rec.seq) in seqs or str(rec.seq.reverse_complement()) in seqs):
			dups += 1; continue
		# check if sequence is low quality
		elif quality_filter(rec, args):
			low_qual += 1; continue
		# keep seq
		else:
			outfile.write('>'+str(read_id)+'\n'+str(rec.seq[0:args['read_length']])+'\n')
			read_id += 1
			if args['filter_dups']: seqs.add(str(sequence))
			if read_id == args['nreads']: break
	# report summary
	if read_id == 0:
		clean_up(paths)
		sys.exit("\tError! No reads remaining after filtering!")
	else:
		args['sampled_reads'] = read_id
	if args['verbose']:
		print ("\t%s reads shorter than %s bp and skipped" % (too_short, args['read_length']))
		print ("\t%s low quality reads found and skipped" % low_qual)
		print ("\t%s duplicate reads found and skipped" % dups)
		print ("\t%s reads sampled from seqfile" % read_id)

def search_seqs(args, paths):
	""" Search high quality reads against marker genes using RAPsearch2 """
	if args['verbose']:
		print ("Searching reads against marker proteins...")
	devnull = open('/dev/null')
	arguments = {'rapsearch': paths['rapsearch'], 'reads':paths['tempfile'], 'db':paths['db'], 'out':paths['tempfile'], 'threads':args['threads']}
	command = "%(rapsearch)s -q %(reads)s -d %(db)s -o %(out)s -z %(threads)s -e 1 -t n -p f -b 0" % arguments
	process = subprocess.Popen(command, shell=True, stdout=devnull, stderr=devnull)
	retcode = process.wait()
	if retcode == 0:
		hits = []
		for line in open(paths['tempfile']+'.m8'):
			if line[0] != '#': hits.append(line.split()[0])
		distinct_hits = len(set(hits))
		if args['verbose']:
			print ("\t%s reads hit marker proteins" % str(distinct_hits))
	else:
		clean_up(paths)
		sys.exit("\tDatabase search has exited with an error!")

def parse_rapsearch(m8):
	""" Yield formatted record from RAPsearch2 m8 file """
	formats = {0:str, 1:str, 2:float, 3:float, 4:float, 5:float, 6:float, 7:float, 8:float, 9:float, 10:float, 11:float}
	fields = {0:'query',1:'target',2:'pid',3:'aln',4:'mis',5:'gaps',6:'qstart',7:'qend',8:'tstart',9:'tend',10:'evalue',11:'score'}
	for line in open(m8):
		if line[0] == '#': continue
		else: yield dict([(fields[index], formats[index](value)) for index, value in enumerate(line.rstrip().split())])

def alignment_coverage(aln_record):
	""" Find the aln cov between query and target
		Do not penalize reads for aligning at gene boundaries
	"""
	r = aln_record
	# convert dna aln coords to prot aln coords for query
	query_len = float(r['query_len'])/3
	query_start_dna, query_stop_dna = sorted([r['qstart'], r['qend']])
	frame = query_start_dna % 3 if query_start_dna % 3 in [1,2] else 3
	query_start = (query_start_dna + 3 - frame)/3
	query_stop = (query_stop_dna + 1 - frame)/3
	# convert aln coords to 1 indexed for target
	target_start, target_stop = sorted([r['tstart'] + 1, r['tend'] + 1])
	# compute aln coverage
	x = min(query_start - 1, target_start - 1)
	y = r['aln']
	z = min(query_len - query_stop, r['target_len'] - target_stop)
	maxaln = x + y + z
	return r['aln'] / maxaln

def alignment_filter(aln_record, optpars):
	""" Determine whether or not ro filter alignment """
	r = aln_record
	if alignment_coverage(aln_record) < optpars[r['target_fam']]['min_cov']:
		return True
	elif r['score'] < optpars[r['target_fam']]['min_score']:
		return True
	elif r['pid'] > optpars[r['target_fam']]['max_aaid']:
		return True
	else:
		return False

def classify_reads(args, paths):
	""" Find best hit for each read and return dictionary {read_id:[fam_id, aln, cov, score]} """
	if args['verbose']:
		print ("Filtering hits...")
	# read in lookups
	optpars  = find_opt_pars(paths['params'], args['read_length'])
	gene2fam = read_dic(paths['fams'], header=False, dtype='char')
	gene2len = read_dic(paths['genelen'], header=False, dtype='float')
	# find best hits
	total_results = 0
	best_hits = {}
	for r in parse_rapsearch(paths['tempfile']+'.m8'):
		total_results += 1
		r['query_len'] = args['read_length']
		r['target_fam'] = gene2fam[r['target']]
		r['target_len'] = gene2len[r['target']]
		if alignment_filter(r, optpars):
			continue
		elif r['query'] not in best_hits:
			best_hits[r['query']] = [r['target_fam'], r['aln'], r['aln']/r['target_len'], r['score']]
		elif best_hits[r['query']][-1] < r['score']:
			best_hits[r['query']] = [r['target_fam'], r['aln'], r['aln']/r['target_len'], r['score']]
	# report summary
	if len(best_hits) == 0:
		clean_up(paths)
		sys.exit("\tError: No hits to marker proteins - cannot estimate genome size! Rerun program with more reads.")
	if args['verbose']:
		print ("\t%s reads assigned to a marker protein" % len(best_hits))
	return best_hits

def aggregate_hits(args, paths, best_hits):
	""" Count hits to each gene family from and return dictionary of results {fam_id:hits} """
	optpars = find_opt_pars(paths['params'], args['read_length'])
	agg_hits = {}
	for fam_id, aln, cov, score in best_hits.values():
		aln_stat = optpars[fam_id]['aln_stat']
		if fam_id not in agg_hits:
			agg_hits[fam_id] = 1.0 if aln_stat == 'hits' else cov if aln_stat == 'cov' else aln
		else:
			agg_hits[fam_id] += 1.0 if aln_stat == 'hits' else cov if aln_stat == 'cov' else aln
	return agg_hits

def estimate_average_genome_size(args, paths, agg_hits):
	"""
		Use linear models to estimate AGS based on the number of hits/bp to each marker (j):
			Rate_j = hits_j/bp
			AGS_j = Coefficient_j/Rate_j
			Model coefficients are listed in <p_coeffs>
		Remove outlier predictions from [AGS_1, AGS_1, ..., AGS_30]
		Take a weighted average across remaining predictions:
			AGS = Sum{j=1 to 30} AGS_j * Weight_j
	"""
	if args['verbose']:
		print ("Computing average genome size...")
	# read in model coefficients and weights
	coeffs  = read_dic(paths['coeffs'],  header=False, dtype='float')
	weights = read_dic(paths['weights'], header=False, dtype='float')
	# predict average genome size for each marker
	estimates = {}
	for fam_id, hits in agg_hits.items():
		coeff = coeffs['_'.join([str(args['read_length']),fam_id])]
		rate = hits/(args['sampled_reads'] * args['read_length'])
		if rate == 0: continue
		else: estimates[fam_id] = coeff/rate
	# take weighted average across predictions & remove outliers   
	est_ags = 0
	mad_estimate = mad(list(estimates.values()))
	median_estimate = median(list(estimates.values()))
	sum_weights = 0
	for fam_id, estimate in estimates.items():
		if abs(estimate - median_estimate) >= mad_estimate:
			continue
		else:
			weight = weights['_'.join([str(args['read_length']),fam_id])]
			est_ags += estimate * weight
			sum_weights += weight
	est_ags = est_ags/sum_weights
	# report results
	if args['verbose']:
		print ("\t%s bp" % str(round(est_ags, 2)))
	return est_ags

def report_results(args, est_ags):
    """ Write estimated average genome size to tab delimited text file """
    outfile = open(args['outfile'], 'w')
    header = ['reads_sampled', 'read_length', 'avg_size']
    outfile.write('\t'.join(header)+'\n')
    record = [args['sampled_reads'], args['read_length'], est_ags]
    outfile.write('\t'.join([str(x) for x in record])+'\n')
    outfile.close()

def clean_up(paths):
	""" Remove all temporary files """
	temp_files = []
	temp_files.append(paths['tempfile'])
	temp_files.append(paths['tempfile']+'.m8')
	temp_files.append(paths['tempfile']+'.aln')
	for file in temp_files:
		if os.path.isfile(file):
			os.remove(file)

def run_pipeline(args):

	# Print copyright
	if 'verbose' in args and args['verbose']: print_copyright()
	
	# Make sure OS is Linux or Darwin
	check_os()

	# Fetch data paths and make sure files exist
	paths = get_relative_paths()
	check_paths(paths)

	# Impute any missing arguments/options, sanity check, print to stdout
	impute_missing_args(args)
	check_arguments(args)
	if args['verbose']: print_parameters(args)
	
	# Sample and QC reads from seqfile
	process_seqfile(args, paths)
	
	# Search reads against marker gene families
	search_seqs(args, paths)
	
	# Classify reads according to optimal parameters
	best_hits = classify_reads(args, paths)
	
	# Count # hits to each gene family
	agg_hits = aggregate_hits(args, paths, best_hits)
		
	# Remove temporary files
	clean_up(paths)

	# Estimate average genome size
	est_ags = estimate_average_genome_size(args, paths, agg_hits)
	return est_ags, args


