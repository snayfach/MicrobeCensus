#######################################################################################
#   LIBRARIES

try:
    import sys
except Exception:
    print 'Module "sys" not installed'; exit()

try:
    import gzip
except Exception:
    print 'Module "gzip" not installed'; exit()

try:
    import random
except Exception:
    print 'Module "random" not installed'; exit()

try:
    import subprocess
except Exception:
    print 'Module "subprocess" not installed'; exit()
    
try:
    import os
except Exception:
    print 'Module "os" not installed'; exit()
    
try:
    import bz2
except Exception:
    print 'Module "bz2" not installed'; exit()     
    
try:
    import Bio.SeqIO
except Exception:
    print 'Module "Bio.SeqIO" not installed'; exit()
    
try:
    import optparse
except Exception:
    print 'Module "optparse" not installed'; exit()

#######################################################################################
#   FUNCTIONS

def mean_stddev(x):
    """ Calculate standard deviation and mean of data in x
    """
    n = float(len(x))
    mean = sum(x)/n
    std  = sqrt(sum([(a - mean)**2 for a in x])/ n)
    return mean, std

def median(x):
    """ Find the median value from unsorted list of numeric values """
    x = sorted(x)
    if len(x) % 2 == 1:
        return x[(len(x)+1)/2-1]
    else:
        lower = x[len(x)/2-1]; upper = x[len(x)/2]
        return (float(lower + upper))/2

def mad(x, const = 1.48):
    """ Calculate the median absolute deviation (MAD) of x
    """
    return const * median([abs(i - median(x)) for i in x])

def clean_up(files):
    """Remove all files in <files>
    """
    for file in files:
        if os.path.isfile(file):
            os.remove(file)

def auto_detect_read_length(p_reads, file_type, p_read_len):
    """ Auto detect median read length from first 10K reads in p_reads
        All reads will be trimmed to this length
    """
    valid_lengths = sorted(read_list(p_read_len, header=True, dtype='int'))
    read_lengths = []
    ext = p_reads.split('.')[-1]
    f_in = gzip.open(p_reads) if ext == 'gz' else bz2.BZ2File(p_reads) if ext == 'bz2' else open(p_reads)
    for index, record in enumerate(Bio.SeqIO.parse(f_in, file_type)):
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
    """ Auto detect file type [fasta or fastq] of p_reads
    """
    ext = p_reads.split('.')[-1]
    f_in = gzip.open(p_reads) if ext == 'gz' else bz2.BZ2File(p_reads) if ext == 'bz2' else open(p_reads) 
    for line in f_in:
        if line[0] == '>': return 'fasta'
        else: return 'fastq'

def auto_detect_fastq_format(p_reads, max_depth):
    """ Auto detect FASTQ file format (sanger, solexa, or illumina)
        For details: http://en.wikipedia.org/wiki/FASTQ_format
    """
    print 'Detecting FASTQ format...'
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
                print '\t', formats.keys()[0]
                return formats.keys()[0]
            elif len(formats) == 0:
                print '\t', 'Unrecognized character in quality string:', line.rstrip()
                exit()
        i = 0; j += 1
        if j == max_depth:
            break
    # guess at format if not detected
    guess = random.sample(formats.keys(), 1)[0]
    print '\t', 'Ambiguous quality encoding, guessing:', guess
    return guess

def process_fasta(p_reads, p_wkdir, nreads, read_length, filter_dups, max_unknown):
    """ Sample <nreads> of <read_length> from <p_reads>
        Write reads to tmpfile in <p_wkdir>
        Return path to tmpfile, sampled read_ids
    """
    print 'Processing sequences...'
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
    for record in Bio.SeqIO.parse(f_in, "fasta"):
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
    print '\t', too_short, 'reads shorter than specified read length and skipped'
    if filter_dups:
        print '\t', dups, 'duplicate reads found and skipped'
    if max_unknown < 1.0:
        print '\t', low_qual, 'low quality reads found and skipped'
    if read_id == 0:
        print '\t', '** Error: No reads remaining after filtering!'
        if keep_tmp is False: clean_up([p_out])
        sys.exit()
    else:
        print '\t', read_id, str(read_length) + 'bp', 'reads sampled from seqfile'
    # return file name, read ids
    return p_out, read_ids

def process_fastq(p_reads, p_wkdir, nreads, read_length, mean_quality, min_quality, filter_dups, max_unknown, fastq_format):
    """ Sample <nreads> of <read_length> from <p_reads> with bases of <min_quality>
        Write reads to tmpfile in <p_wkdir>
        Return path to tmpfile, sampled read_ids
    """
    print 'Processing sequences...'
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
    for record in Bio.SeqIO.parse(f_in, encoding):
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
    print '\t', too_short, 'reads shorter than specified read length and skipped'
    if filter_dups:
        print '\t', dups, 'duplicate reads found and skipped'
    if min_quality > -5 or mean_quality > -5 or max_unknown < 1.0:
        print '\t', low_qual, 'low quality reads found and skipped'
    if read_id == 0:
        print '\t', '** Error: No reads remaining after filtering!'
        if keep_tmp is False: clean_up([p_out])
        sys.exit()
    else:
        print '\t', read_id, str(read_length) + 'bp', 'reads sampled from seqfile'
    # return file name, read ids
    return p_out, read_ids

def search_seqs(reads, db, rapsearch, threads, keep_tmp, p_params):
    """ Searches fasta file <reads> against <db> using <threads> 
        Writes search results to temporary file in dir of <reads>
        Returns path to <reads>.m8 and <reads>.aln
    """
    print 'Searching marker proteins...'
    out = reads
    devnull = open('/dev/null')
    arguments = {'rapsearch': rapsearch, 'reads':reads, 'db':db, 'out':out, 'threads':threads}
    command = "%(rapsearch)s -q %(reads)s -d %(db)s -o %(out)s -z %(threads)s -e 1 -t n -p f -b 0" % arguments
    process = subprocess.Popen(command, shell=True, stdout=devnull, stderr=devnull)
    retcode = process.wait()
    if retcode == 0:
        if keep_tmp is False: clean_up([reads])
        return (out+'.m8', out+'.aln')
    else:
        print '** Error: Database search has exited with an error!'
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
    """ Parse and filter results, classify reads into families
        Use class parameters specific to read length and marker family
        Returns dic {read_id:[fam_id, aln, cov, score]}
    """
    print 'Filtering alignments...'
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
        print '** Error: No hits to marker proteins - cannot estimate genome size! Rerun program with more reads.'
        sys.exit()
    else:
        print '\t', str(len(best_hits))+'/'+str(n_reads_sampled), 'reads classified into a single-copy gene family'
    return best_hits

def resample_hits(best_hits, read_ids):
    """ Find best hit for each bootstrapped query
        Returns list [fam_id, aln, cov, score]
    """
    # randomly sample read_ids for bootstrapping
    random_reads = {}
    for i in range(len(read_ids)):
        read_id = random.sample(read_ids, 1)[0]
        if read_id not in random_reads:
            random_reads[read_id] = 1
        else:
            random_reads[read_id] += 1
    # find best hits
    resampled_hits = {}
    for rand_read_id, n in random_reads.iteritems():
        if rand_read_id in best_hits:
            for i in range(n):
                resampled_read_id = rand_read_id + '_' + str(n)
                resampled_hits[resampled_read_id] = best_hits[rand_read_id]
        else:
            continue
    # check that at least 1 filtered read was selected
    if len(resampled_hits) == 0:
        print '**Error: No hits to marker proteins - cannot estimate genome size! Rerun program with more reads.'
        sys.exit()
    return resampled_hits
           
def aggregate_hits(best_hits, read_length, p_params):
    """ Get sum of alignments to marker families from best hits 
        Sum the alignment statistic (hits, cov, aln) specific to read length and marker family
        Returns agg_hits[fam_id] = hits
    """
    optpars = find_opt_pars(p_params, read_length)
    agg_hits = {}
    for fam_id, aln, cov, score in best_hits.values():
        aln_stat = optpars[fam_id][-1]
        if fam_id not in agg_hits: agg_hits[fam_id] = 1.0 if aln_stat == 'hits' else cov if aln_stat == 'cov' else aln
        else: agg_hits[fam_id] += 1.0 if aln_stat == 'hits' else cov if aln_stat == 'cov' else aln
    return agg_hits

def pred_genome_size(agg_hits, n_reads_sampled, read_length, p_coeffs, p_weights):
    """ Return the average genome size of input metagenome
    """
    print 'Computing average genome size..'
    # read in model coefficients and weights
    coeffs  = read_dic(p_coeffs,  header=False, dtype='float')
    weights = read_dic(p_weights, header=False, dtype='float')
    
    # predict genome sizes for each marker
    preds = {}
    for fam_id in agg_hits.keys():
        coeff   = coeffs[str(read_length) + '_' + fam_id]
        hits    = agg_hits[fam_id]
        rate    = hits/(n_reads_sampled * read_length)
        if rate == 0:
            print '\t**Warning: no hits to gene', fam_id + '. Skipping.'
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

def bootstrapping(best_hits, read_ids, read_length, n_reads_sampled, nboot):
    """ 
        For each boostrap iteration:
            1. Sample reads with replacement from <best_hits> 
            2. Aggregate hits to each family
            3. Estimate coverage and genome size using each family
            4. Take a weighted average of estimates
        Return the mean coverage and genome size across bootstrap iterations
    """
    print 'Bootstrapping...'
    boot_coverage = []
    boot_avg_size = []
    for n in range(nboot):
        resampled_hits                 = resample_hits(best_hits, read_ids)
        agg_hits                       = aggregate_hits(best_hits, read_length)
        avg_size                       = pred_genome_size(agg_hits, n_reads_sampled, read_length)
        boot_avg_size.append(avg_size)
    mean_avg_size, stddev_avg_size     = mean_stddev(boot_avg_size)
    return mean_avg_size, stddev_avg_size

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