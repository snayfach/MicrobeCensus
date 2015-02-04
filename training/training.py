# MicrobeCensus - estimation of average genome size from shotgun sequence data
# Copyright (C) 2013-2014 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

import subprocess, gzip, Bio.SeqIO, os, numpy
from time import sleep
from multiprocessing import Process, Queue

def parallel_subprocess(command, args_list, threads):
	""" Run command using subprocess in parallel """
	processes = []
	for arguments in args_list: # run command for each set of args in args_list
		p = subprocess.Popen(command % arguments, shell=True, stdout=open('/dev/null', 'w'), stderr=subprocess.STDOUT)
		processes.append(p)
		while len(processes) >= threads: # control number of active processes
			sleep(1)
			indexes = []
			for index, process in enumerate(processes):
				if process.poll() is None: # has process return retcode?
					indexes.append(index)  # keep process
			processes = [processes[i] for i in indexes] # update list of processes
	while None in [p.poll() for p in processes]: # wait until there are no active processes
		sleep(1)

def parallel_function(function, args_list, threads):
	""" Run function using multiple threads """
	processes = []
	for arguments in args_list: # run function for each set of args in args_list
		p = Process(target=function, kwargs=arguments)
		processes.append(p)
		processes[-1].start()
		while len(processes) >= threads: # control number of active processes
			sleep(1)
			indexes = []
			for index, process in enumerate(processes):
				if process.is_alive(): # keep processes that are still alive
					indexes.append(index)
			processes = [processes[i] for i in indexes] # update list of processes
	while any([p.is_alive() for p in processes]): # wait until there are no active processes
		sleep(1)

def parallel_return_function(function, args_list, threads):
	""" Run function using multiple threads """
	values = []
	processes = []
	queue = Queue()
	for arguments in args_list: # run function for each set of args in args_list
		arguments['queue'] = queue # append queue to list of args
		p = Process(target=function, kwargs=arguments)
		processes.append(p)
		processes[-1].start()
		while len(processes) >= threads: # control number of active processes
			indexes = []
			for index, process in enumerate(processes):
				if process.is_alive(): # keep processes that are still alive
					indexes.append(index)
				else: # return values from processes that are finished
					values.append(queue.get())
			processes = [processes[i] for i in indexes] # update list of processes
	# wait until there are no active processes
	while len(processes) > 0:
		indexes = []
		for index, process in enumerate(processes):
			if process.is_alive(): # keep processes that are still alive
				indexes.append(index)
			else: # return values from processes that are finished
				values.append(queue.get())
		processes = [processes[i] for i in indexes] # update list of processes
	return values

def xvalidation(pars, xfolds, genome_names, rates, genome2size, queue):
	""" Perform xfold cross validation to estimate test error for each set of parameters """
	error = []
	# loop over folds
	for i in range(1, xfolds + 1):
		# split data into training and testing
		train_i, test_i = xfold_indexes(len(genome_names), xfolds, i)
		train_genomes   = [genome_names[i] for i in train_i]
		train_rates     = [rates[i] for i in train_i]
		test_genomes    = [genome_names[i] for i in test_i]
		test_rates      = [rates[i] for i in test_i]
		# train model
		prop_constant = estimate_proportionality_constant(train_genomes, train_rates, genome2size)
		# test model
		error += test_error(train_genomes, train_rates, prop_constant, genome2size)
	queue.put([pars, numpy.median(error)])

def genome_sizes(genomes_dir):
	""" Return sizes (bp) of all genomes in directory """
	genome2size = {}
	for genome_file in os.listdir(genomes_dir):
		if genome_file[-7:] != '.fna.gz': sys.exit('Genomes must have .fna.gz extension')
		genome_path = os.path.join(genomes_dir, genome_file)
		genome_name = genome_file[:-7]
		genome2size[genome_name] = compute_seq_len(genome_path)
	return genome2size

def library_sizes(reads_dir):
	""" Return sizes (bp) of all libraries in directory """
	library2size = {}
	for read_length in os.listdir(reads_dir):
		for reads_file in os.listdir(os.path.join(reads_dir, read_length)):
			if reads_file[-9:] != '-reads.fa': continue
			reads_path = os.path.join(os.path.join(reads_dir, read_length), reads_file)
			genome_name = reads_file[:-9]
			library2size[(read_length, genome_name)] = compute_seq_len(reads_path)
	return library2size

def xfold_indexes(n, x, i):
	""" Split n data points into training and testing groups
		
		n = total data points
	    x = folds for cross-validation (ex: 10-fold)
		i = iteration (1, 2, ..., x)
	"""
	fold_size = n/x
	test_start = fold_size * i - fold_size
	test_stop = test_start + fold_size - 1
	test_indexes = range(test_start, test_stop + 1)
	train_indexes = range(n)
	for j in test_indexes:
		train_indexes.remove(j)
	return train_indexes, test_indexes

def estimate_proportionality_constant(genome_names, rates, genome2size):
	""" Estimate proportionality constant for AGS estimation
		S_exp = C_est/R_obs -> C_est = S_exp * R_obs
			S_exp: expected genome size
			C_est: estimated proportionality constant
			R_obs: observed classification rate
		Return the median proportionality constant across simulations
	"""
	constants = []
	for genome_name, rate in zip(genome_names, rates):
		genome_size = genome2size[genome_name]
		constants.append(genome_size * rate)
	return numpy.median(constants)

def test_error(genome_names, rates, prop_constant, genome2size):
	""" Compute prediction error for test genomes using estimated proportionality constant
		S_est = C_est/R_obs
			S_est: estimated genome size, 
			C_est: estimated proportionality constant
			R_obs: observed classification rate
		E = 100 * |S_est - S_exp|/S_exp
			E: percent error
			S_est: estimated genome size
			S_exp: expected genome size
		Return E for all test_genomes
	"""
	error = []
	for genome_name, rate in zip(genome_names, rates):
		if rate == 0: # if rate == 0, genome size cannot be estimated
			error.append(999)
		else:
			pred_size = prop_constant/rate
			true_size = genome2size[genome_name]
			error.append(100 * abs(pred_size - true_size)/true_size)
	return error

def find_opt_pars(xval_error):
	opt_pars = {}
	for pars, error in xval_error:
		read_length, fam, min_score, max_pid, aln_cov, rate_type = pars
		if (read_length, fam) not in opt_pars:
			opt_pars[read_length, fam] = {'pars': (min_score, max_pid, aln_cov, rate_type), 'error':error}
		elif error < opt_pars[read_length, fam]['error']:
			opt_pars[read_length, fam] = {'pars': (min_score, max_pid, aln_cov, rate_type), 'error':error}
	return opt_pars

def store_rates(hits_dir, library2size):
	""" Store classification rates across to each marker across libraries """
	rates = {}
	# loop over hits directories
	for read_length in os.listdir(hits_dir):
		rates[read_length] = {}
		for hits_file in os.listdir(os.path.join(hits_dir, read_length)):
			if hits_file[-5:] != '.hits': continue
			hits_path = os.path.join(os.path.join(hits_dir, read_length), hits_file)
			genome_name = hits_file[0:-5]
			for r in parse_hits(hits_path):
				fam = r['fam']
				map_pars = (r['min_score'], r['max_pid'], r['aln_cov'])
				# initialize keys
				if fam not in rates[read_length]:
					rates[read_length][fam] = {}
				if map_pars not in rates[read_length][fam]:
					rates[read_length][fam][map_pars] = {}
					rates[read_length][fam][map_pars]['genome_names'] = []
					rates[read_length][fam][map_pars]['rate_hits'] = []
					rates[read_length][fam][map_pars]['rate_aln'] = []
					rates[read_length][fam][map_pars]['rate_cov'] = []
				# store data
				rates[read_length][fam][map_pars]['genome_names'].append(genome_name)
				rates[read_length][fam][map_pars]['rate_hits'].append(r['count_hits']/library2size[(read_length, genome_name)])
				rates[read_length][fam][map_pars]['rate_aln'].append(r['count_aln']/library2size[(read_length, genome_name)])
				rates[read_length][fam][map_pars]['rate_cov'].append(r['count_cov']/library2size[(read_length, genome_name)])
	return rates

def parse_hits(p_in):
	""" Yield data-type formatted dictionary of records from .hits file """
	formats = {0:str,  1:float,    2:float,    3:float,      4:float,       5:float,      6:float      }
	fields  = {0:'fam',1:'aln_cov',2:'max_pid',3:'min_score',4:'count_hits',5:'count_aln',6:'count_cov'}
	f_in = open(p_in)
	next(f_in)
	for line in f_in:
		x = line.rstrip().split()
		yield dict( [ (fields[index], formats[index](value)) for index, value in enumerate(x) ] )

def parse_rapsearch(p_in):
	""" Yield data-type formatted dictionary of records from RAPsearch2 m8 file """
	formats = {0:str, 1:str, 2:float, 3:int, 4:float, 5:float, 6:float, 7:float, 8:float, 9:float, 10:float, 11:float}
	fields = {0:'query',1:'target',2:'pid',3:'aln',4:'mismatches',5:'gaps',6:'qstart',7:'qend',8:'tstart',9:'tend',10:'log_evalue',11:'score'}
	f_in = open(p_in)
	for line in f_in:
		if line[0] == '#': continue
		x = line.rstrip().split()
		z = dict( [ (fields[index], formats[index](value)) for index, value in enumerate(x) ] )
		yield z

def compute_seq_len(p_in):
	""" Return total sequence length of input fasta file """
	size = 0
	f_in = gzip.open(p_in) if p_in[-3:] == '.gz' else open(p_in)
	for r in Bio.SeqIO.parse(f_in, 'fasta'):
		size += len(r.seq)
	return size

def read_hits(p_in, gene2fam):
	""" Read in hits into list """
	my_hits = []
	for r in parse_rapsearch(p_in):
		# convert query start/stop to amino acid space
		query_start_dna, query_stop_dna = sorted([r['qstart'], r['qend']])
		frame = query_start_dna % 3 if query_start_dna % 3 in [1,2] else 3
		query_start = (query_start_dna + 3 - frame)/3
		query_stop = (query_stop_dna + 1 - frame)/3
		# make target alignment coords 1 indexed instead of 0 indexed
		target_start, target_stop = sorted([r['tstart'] + 1, r['tend'] + 1])
		target_fam = gene2fam[r['target']]
		# store results
		hit = [ r['query'], r['target'], target_fam, r['pid'], r['aln'], query_start, query_stop, target_start, target_stop, r['score'] ]
		my_hits.append(hit)
	return my_hits

def aln_filter(hits, aln_cov, read_length, gene2len):
	""" Filter hits by alignment coverage """
	hits_aln_filt = []
	for hit in hits:
		# parse record
		read_id, target_gene, target_fam, pid, aln, query_start, query_stop, target_start, target_stop, score = hit
		query_len = float(read_length)/3
		target_len = gene2len[target_gene]
		# determine maximum aln length given gene boundary
		x = min(query_start - 1, target_start - 1)
		y = aln
		z = min(query_len - query_stop, target_len - target_stop)
		maxaln = x + y + z
		# skip records with insufficient alignment coverage
		if aln/maxaln < aln_cov:
			continue
		# store results
		hits_aln_filt.append(hit)
	return hits_aln_filt

def pid_filter(hits, max_pid):
	""" Filter hits by percent identity """
	hits_pid_filt = []
	for hit in hits:
		read_id, target_gene, target_fam, pid, aln, query_start, query_stop, target_start, target_stop, score = hit
		if pid > max_pid:
			continue
		else:
			hits_pid_filt.append(hit)
	return hits_pid_filt

def score_filter(hits, min_score):
	""" Filter hits by bit-score """
	hits_score_filt = []
	for hit in hits:
		read_id, target_gene, target_fam, pid, aln, query_start, query_stop, target_start, target_stop, score = hit
		if score < min_score:
			continue
		else:
			hits_score_filt.append(hit)
	return hits_score_filt

def find_best_hits(hits):
	""" Find top-scoring hit for each read """
	best_hits = {}
	for hit in hits:
		read_id, target_gene, target_fam, pid, aln, query_start, query_stop, target_start, target_stop, score = hit
		if read_id not in best_hits:
			best_hits[read_id] = hit
		elif best_hits[read_id][-1] < score:
			best_hits[read_id] = hit
	return best_hits.values()

def aggregate_hits(hits, fams, gene2len):
	""" Aggregate hits to each gene family """
	fam_2_hits = {}
	for fam in fams:
		fam_2_hits[fam] = {'hits':0, 'cov':0, 'aln':0}
	for hit in hits:
		read_id, target_gene, target_fam, pid, aln, query_start, query_stop, target_start, target_stop, score = hit
		fam_2_hits[target_fam]['hits'] += 1
		fam_2_hits[target_fam]['cov']  += float(aln)/gene2len[target_gene]
		fam_2_hits[target_fam]['aln']  += aln
	return fam_2_hits

def classify_reads(p_in, p_out, aln_covs, max_pids, min_scores, gene2len, gene2fam, fams, read_length):
	""" Perform grid search, classifying reads across all combinations of maping parameters (aln_covs, max_pids, min_scores)
	"""
	# open output
	f_out = open(p_out, 'w')
	f_out.write('\t'.join(['fam', 'aln_cov', 'max_pid', 'min_score', 'count_hits', 'count_aln', 'count_cov'])+'\n')
	# read in hits
	hits = read_hits(p_in, gene2fam)
	# alignment coverage filter
	for aln_cov in aln_covs:
		hits_aln_filt = aln_filter(hits, aln_cov, read_length, gene2len)
		# percent identity filter
		for max_pid in max_pids:
			hits_pid_filt = pid_filter(hits_aln_filt, max_pid)
			# score filter
			for min_score in min_scores:
				hits_score_filt = score_filter(hits_pid_filt, min_score)
				# find best hits
				best_hits = find_best_hits(hits_score_filt)
				# count hits to each marker familes
				for fam, stats in aggregate_hits(best_hits, fams, gene2len).iteritems():
					record = [str(x) for x in [fam, aln_cov, max_pid, min_score, stats['hits'], stats['aln'], stats['cov']]]
					# write results to p_out
					f_out.write('\t'.join(record)+'\n')

def drange(start, stop, step):
	""" Decimal range """
	my_range = []
	r = start
	while r < stop:
		my_range.append(r)
		r += step
	return my_range