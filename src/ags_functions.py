# MicrobeCensus - estimation of average genome size from shotgun sequence data
# Copyright (C) 2013-2014 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

#######################################################################################
#   LIBRARIES

try:
    import sys
except Exception:
    print ('Module "sys" not installed'); exit()

try:
    import gzip
except Exception:
    print ('Module "gzip" not installed'); exit()

try:
    import random
except Exception:
    print ('Module "random" not installed'); exit()

try:
    import subprocess
except Exception:
    print ('Module "subprocess" not installed'); exit()
    
try:
    import os
except Exception:
    print ('Module "os" not installed'); exit()
    
try:
    import bz2
except Exception:
    print ('Module "bz2" not installed'); exit()
    
try:
    from Bio.SeqIO import parse
except Exception:
    print ('Module "Bio.SeqIO" not installed'); exit()
    
try:
    import optparse
except Exception:
    print ('Module "optparse" not installed'); exit()

try:
	from numpy import median
	from numpy import mean
	from numpy import std
except Exception:
    print ('Module "numpy" not installed'); exit()


#######################################################################################
#   FUNCTIONS

def mad(x, const = 1.48):
    """ Calculate the median absolute deviation (MAD) of <x>
    """
    return const * median([abs(i - median(x)) for i in x])

def clean_up(files):
    """Remove all files in <files>
    """
    for file in files:
        if os.path.isfile(file):
            os.remove(file)

def auto_detect_read_length(p_reads, file_type, p_read_len):
    """ Find median read length from first 10K reads in <p_reads>
		Value rounded down to nearest legal read length listed in <p_read_len>
    """
    valid_lengths = sorted(read_list(p_read_len, header=False, dtype='int'))
    read_lengths = []
    ext = p_reads.split('.')[-1]
    f_in = gzip.open(p_reads) if ext == 'gz' else bz2.BZ2File(p_reads) if ext == 'bz2' else open(p_reads)
    for index, record in enumerate(parse(f_in, file_type)):
        if index == 10000: break
        read_lengths.append(len(record.seq))
    median_read_length = int(median(read_lengths))
    if median_read_length < valid_lengths[0]:
        sys.exit('Median read length is %s. Cannot compute AGS using reads shorter than 50 bp.' % str(median_read_length))
    for index, read_length in enumerate(valid_lengths):
        if read_length > median_read_length:
            return valid_lengths[index-1]
    return valid_lengths[-1]

def auto_detect_file_type(p_reads):
    """ Detect file type [fasta or fastq] of <p_reads>
    """
    ext = p_reads.split('.')[-1]
    f_in = gzip.open(p_reads) if ext == 'gz' else bz2.BZ2File(p_reads) if ext == 'bz2' else open(p_reads) 
    for line in f_in:
        if line[0] == '>': return 'fasta'
        else: return 'fastq'

def auto_detect_fastq_format(p_reads):
	""" Detect quality score encoding of <p_reads> file (sanger, solexa, or illumina)
		For details: http://en.wikipedia.org/wiki/FASTQ_format
	"""
	max_depth = 1000000
	# format specific characters
	sanger   = set(list("""!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"""))
	solexa   = set(list(""";<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"""))
	illumina = set(list("""@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefgh"""))
	formats  = {'sanger': sanger, 'solexa': solexa, 'illumina': illumina}
	# open input stream
	ext = p_reads.split('.')[-1]
	f_in = gzip.open(p_reads) if ext == 'gz' else bz2.BZ2File(p_reads) if ext == 'bz2' else open(p_reads)
	# loop over quality scores
	i = 0; j = 0
	for line in f_in:
		i += 1
		if i != 4: continue
		for q in line.rstrip():
			# look for incompatible formats
			for format in formats.keys():
				if q not in formats[format]:
					del formats[format]
			if len(formats) == 1:
				return formats.keys()[0]
			elif len(formats) == 0:
				print ('\t', 'Unrecognized character in quality string:', line.rstrip())
				exit()
		i = 0; j += 1
		if j == max_depth:
			break
	# guess at format if not detected
	guess = random.sample(formats.keys(), 1)[0]
	return guess

def process_fasta(p_reads, p_wkdir, nreads, read_length, filter_dups, max_unknown, keep_tmp):
    """ Sample <nreads> of <read_length> from <p_reads>
        Write reads to tmpfile in <p_wkdir>
        Return path to tmpfile, sampled read_ids
    """
    print ('Sampling & trimming reads...')
    # open input, output files
    ext = p_reads.split('.')[-1]
    f_in = gzip.open(p_reads) if ext == 'gz' else bz2.BZ2File(p_reads) if ext == 'bz2' else open(p_reads)
    p_out = os.path.join(p_wkdir, os.path.basename(p_reads).rstrip('.gz')+'.tmp')
    f_out = open(p_out, 'w')
    # loop over sequences
    read_ids = []
    read_id  = 0
    dups = 0
    too_short = 0
    low_qual = 0
    seqs = set([])
    for record in parse(f_in, "fasta"):
        # parse record
        id = str(record.id)
        sequence = record.seq
        my_seq = sequence[0:read_length]
        # record sequence if enough high quality bases remain
        if len(sequence) < read_length:
            too_short += 1
            continue
        # check if sequence is a duplicate
        elif filter_dups is True and (str(sequence) in seqs or str(sequence.reverse_complement()) in seqs):
            dups += 1
            continue
        # check proportion of ambiguous base calls (N)
        elif sum([1 if base == 'N' else 0 for base in list(my_seq)])/float(read_length) > max_unknown:
            low_qual += 1
            continue
        # keep seq
        else:
            f_out.write('>' + str(read_id) + '\n' + str(sequence[0:read_length]) + '\n')
            read_ids.append(read_id)
            read_id += 1
            if filter_dups is True:
                seqs.add(str(sequence))
        # stop when nreads written
        if read_id == nreads:
            break
    # print status: # of reads sampled, filtered
    print ('\t%s reads shorter than %s bp and skipped' % (too_short, read_length))
    if filter_dups:
        print ('\t%s duplicate reads found and skipped' % dups)
    if max_unknown < 1.0:
        print ('\t%s low quality reads found and skipped' % low_qual)
    if read_id == 0:
        print ('\t** Error: No reads remaining after filtering!')
        if keep_tmp is False: clean_up([p_out])
        sys.exit()
    else:
        print ('\t%s reads sampled from seqfile' % read_id)
    # return file name, read ids
    return p_out, read_ids

def process_fastq(p_reads, p_wkdir, nreads, read_length, mean_quality, min_quality, filter_dups, max_unknown, fastq_format, keep_tmp):
    """ Sample <nreads> of <read_length> from <p_reads> with bases of <min_quality>
        Write reads to tmpfile in <p_wkdir>
        Return path to tmpfile, sampled read_ids
    """
    print ('Sampling & trimming reads...')
    # open input, output files
    ext = p_reads.split('.')[-1]
    f_in = gzip.open(p_reads) if ext == 'gz' else bz2.BZ2File(p_reads) if ext == 'bz2' else open(p_reads)
    p_out = os.path.join(p_wkdir, os.path.basename(p_reads).rstrip('.gz')+'.tmp')
    f_out = open(p_out, 'w')
    # determine format
    encoding = "fastq-sanger" if fastq_format == "sanger" else "fastq-solexa" if fastq_format == "solexa" else "fastq-illumina"
    quality = "phred_quality" if fastq_format in ["sanger", "illumina"] else "solexa_quality"
    # loop over sequences
    read_ids = []
    read_id  = 0
    seqs = set([])
    dups = 0
    low_qual = 0
    too_short = 0
    for record in parse(f_in, encoding):
        # parse record
        id = str(record.id)
        sequence = record.seq
        my_qscores = record.letter_annotations[quality][0:read_length]
        my_seq = sequence[0:read_length]
        # check if sequence is long enough
        if len(my_seq) < read_length:
            too_short += 1
            continue        
        # check if sequence is a duplicate
        elif filter_dups is True and (str(sequence) in seqs or str(sequence.reverse_complement()) in seqs):
            dups += 1
            continue
        # check read-level sequence quality
        elif sum(my_qscores)/float(read_length) < mean_quality:
            low_qual += 1
            continue        
        # check base-level sequence quality
        elif any([q < min_quality for q in list(my_qscores)]):
            low_qual += 1
            continue
        # check proportion of ambiguous base calls (N)
        elif sum([1 if base == 'N' else 0 for base in list(my_seq)])/float(read_length) > max_unknown:
            low_qual += 1
            continue
        # keep seq
        else:
            f_out.write('>' + str(read_id) + '\n' + str(my_seq) + '\n')
            read_ids.append(read_id)
            read_id += 1
            if filter_dups is True:
                seqs.add(str(sequence))
        # stop when nreads written
        if read_id == nreads:
            break
    # print status: # of reads sampled, filtered
    print ('\t%s reads shorter than %s bp and skipped' % (too_short, read_length))
    if filter_dups:
        print ('\t%s duplicate reads found and skipped' % dups)
    if min_quality > -5 or mean_quality > -5 or max_unknown < 1.0:
        print ('\t%s low quality reads found and skipped' % low_qual)
    if read_id == 0:
        print ('\t** Error: No reads remaining after filtering!')
        if keep_tmp is False: clean_up([p_out])
        sys.exit()
    else:
        print ('\t%s reads sampled from seqfile' % read_id)
    # return file name, read ids
    return p_out, read_ids

def search_seqs(reads, db, rapsearch, threads, keep_tmp, p_params):
	""" Searches fasta file <reads> against <db> using <threads>
		Writes search results to temporary file in dir of <reads>
		Returns path to <reads>.m8 and <reads>.aln
	"""
	print ('Searching reads against marker proteins...')
	out = reads
	devnull = open('/dev/null')
	arguments = {'rapsearch': rapsearch, 'reads':reads, 'db':db, 'out':out, 'threads':threads}
	command = "%(rapsearch)s -q %(reads)s -d %(db)s -o %(out)s -z %(threads)s -e 1 -t n -p f -b 0" % arguments
	process = subprocess.Popen(command, shell=True, stdout=devnull, stderr=devnull)
	retcode = process.wait()
	if retcode == 0:
		if keep_tmp is False:
			clean_up([reads])
		hits = []
		for line in open(out+'.m8'):
			if line[0] != '#': hits.append(line.split()[0])
		distinct_hits = len(set(hits))
		print ('\t%s reads hit marker proteins' % str(distinct_hits))
		return (out+'.m8', out+'.aln')
	else:
		print ('** Error: Database search has exited with an error!')
		if keep_tmp is False: clean_up([reads, out+'.m8', out+'.aln'])
		sys.exit()

def gmaxaln_cov(aln, read_length, target_len, query_start, query_stop, target_start, target_stop):
    """ Find the aln cov between query and target
        Computes cov based on max possible aln between seqs
        Does not penalize reads for aligning at gene boundaries
    """
    # convert dna aln coords to prot aln coords for query
    query_len = float(read_length)/3
    query_start_dna, query_stop_dna = sorted([query_start, query_stop])
    frame = query_start_dna % 3 if query_start_dna % 3 in [1,2] else 3
    query_start = (query_start_dna + 3 - frame)/3
    query_stop = (query_stop_dna + 1 - frame)/3
    # convert aln coords to 1 indexed for target
    target_start, target_stop = sorted([target_start + 1, target_stop + 1])
    # compute aln coverage
    x = min(query_start - 1, target_start - 1)
    y = aln
    z = min(query_len - query_stop, target_len - target_stop)
    maxaln = x + y + z
    return aln/maxaln

def classify_reads(results, alignments, read_length, p_params, p_gene2fam, p_gene2len, keep_tmp, n_reads_sampled):
    """ Find best hit for each read from <results>
        Use thresholds listed in <p_params> to classify reads into marker families. These are specific to each family and read length
        Return dictionary of classified reads {read_id:[fam_id, aln, cov, score]}
    """
    print ('Filtering hits...')
    # read in lookups
    optpars  = find_opt_pars(p_params, read_length)
    gene2fam = read_dic(p_gene2fam, header=False, dtype='char')
    gene2len = read_dic(p_gene2len, header=False, dtype='float')
    # find best hits
    total_results = 0
    best_hits = {}
    f_in = open(results)
    for line in f_in:
        # parse record
        if line[0] == '#': continue
        total_results += 1
        query = line.rstrip().split()[0]
        target_gene  = line.rstrip().split()[1]
        target_fam   = gene2fam[target_gene]
        target_len   = gene2len[target_gene]
        aaid, aln, mismatches, gaps, query_start, query_stop, target_start, target_stop, eVal, score = [ float(x) for x in line.rstrip().split()[2:] ]
        cov = gmaxaln_cov(aln, read_length, target_len, query_start, query_stop, target_start, target_stop)
        min_cov, max_aaid, min_score, aln_stat = optpars[target_fam]
        target_cov = aln/target_len
        # filter record
        if cov < min_cov or score < min_score or aaid > max_aaid:
            continue
        # store record
        if query not in best_hits:
            best_hits[query] = [target_fam, aln, target_cov, score]
        elif best_hits[query][-1] < score:
            best_hits[query] = [target_fam, aln, target_cov, score]
    # check that at least 1 read passed filters
    if keep_tmp is False: clean_up([results, alignments])
    if len(best_hits) == 0:
        print ('** Error: No hits to marker proteins - cannot estimate genome size! Rerun program with more reads.')
        sys.exit()
    else:
        print ('\t%s reads assigned to a marker protein' % len(best_hits))
    return best_hits
           
def aggregate_hits(best_hits, read_length, p_params):
    """ Sum # of alignments to each marker family from <best_hits>
        Return dictionary of aggregated hits {fam_id:hits}
    """
    optpars = find_opt_pars(p_params, read_length)
    agg_hits = {}
    for fam_id, aln, cov, score in best_hits.values():
        aln_stat = optpars[fam_id][-1]
        if fam_id not in agg_hits: agg_hits[fam_id] = 1.0 if aln_stat == 'hits' else cov if aln_stat == 'cov' else aln
        else: agg_hits[fam_id] += 1.0 if aln_stat == 'hits' else cov if aln_stat == 'cov' else aln
    return agg_hits

def pred_genome_size(agg_hits, n_reads_sampled, read_length, p_coeffs, p_weights):
    """ 
		Use linear models to estimate AGS based on the number of hits/bp to each marker (j):
			Rate_j = hits_j/bp
			AGS_j = Coefficient_j/Rate_j
			Model coefficients are listed in <p_coeffs>
		Remove outlier predictions from [AGS_1, AGS_1, ..., AGS_30]
		Take a weighted average across remaining predictions:
			AGS = Sum{j=1 to 30} AGS_j * Weight_j
			Model weights are listed in <p_weights>
    """
    print ('Computing average genome size..')
    # read in model coefficients and weights
    coeffs  = read_dic(p_coeffs,  header=False, dtype='float')
    weights = read_dic(p_weights, header=False, dtype='float')
    # predict average genome size for each marker
    preds = {}
    for fam_id in agg_hits.keys():
        coeff   = coeffs[str(read_length) + '_' + fam_id]
        hits    = agg_hits[fam_id]
        rate    = hits/(n_reads_sampled * read_length)
        if rate == 0:
            print ('\t**Warning: no hits to gene', fam_id + '. Skipping.')
            continue
        else:
            pred = coeff/rate
            preds[fam_id] = pred
    # take weighted average across predictions & remove outliers   
    pred_size   = 0
    mad_pred    = mad(preds.values())
    median_pred = median(preds.values())
    sum_weights = 0
    for fam_id in preds.keys():
        if abs(preds[fam_id] - median_pred) >= mad_pred:
            continue
        else:
            pred_size   += preds[fam_id] * weights[str(read_length)+'_'+fam_id]
            sum_weights += weights[str(read_length)+'_'+fam_id]
    pred_size = pred_size/sum_weights
    return pred_size

def write_results(p_out, n_reads_sampled, read_length, avg_size, keep_tmp, p_results, p_aln):
    """ Write avg genome size to tab delimited text file
        Clean up remaining temporary files
    """
    # write results
    f_out = open(p_out, 'w')
    header = ['reads_sampled', 'read_length', 'avg_size']
    f_out.write('\t'.join(header)+'\n')
    record = [n_reads_sampled, read_length, avg_size]
    f_out.write('\t'.join([str(x) for x in record])+'\n')
    f_out.close()
    # cleanup temporary files
    if keep_tmp is False: clean_up([p_results, p_aln])

def find_opt_pars(p_optpars, read_length):
    """ Read in optimal parameters for each family at given read length
        Returns optpars[fam_id] = [min_cov, max_aaid, min_score, aln_stat]
    """
    optpars = {}
    with open(p_optpars) as f_in:
        next(f_in)
        for line in f_in:
            fam_id, read_len, min_cov, max_aaid, min_score, aln_stat = line.rstrip().split()
            if read_len != str(read_length): continue
            else: optpars[fam_id] = [float(min_cov), float(max_aaid), float(min_score), aln_stat]
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
    """ Read in simple set from file 
    """
    my_list = []
    f_in = gzip.open(file) if file.split('.')[-1] == 'gz' and file.split('.') > 1 else open(file)
    if header is True: next(f_in)
    for line in f_in:
        value = line.rstrip()
        my_list.append(float(value)) if dtype == 'float' else my_list.append(int(value)) if dtype == 'int' else my_list.append(value)
    return my_list