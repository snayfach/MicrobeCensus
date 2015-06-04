# MicrobeCensus - estimation of average genome size from shotgun sequence data
# Copyright (C) 2013-2015 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

# import libraries
try:
    import random
except Exception:
    print 'Module "random" not installed'; exit()

try:
    import gzip
except Exception:
    print 'Module "gzip" not installed'; exit()

try:
    import os
except Exception:
    print 'Module "os" not installed'; exit()

try:
    import Bio.SeqIO
except Exception:
    print 'Module "Bio.SeqIO" not installed'; exit()
    
try:
    import math
except Exception:
    print 'Module "math" not installed'; exit()


#   define functions
def read_in_genome(p_genome):
	ext = p_genome.split('.')[-1]
	handle = gzip.open(p_genome) if ext == 'gz' else open(p_genome)
	return dict([(r.id, r.seq) for r in Bio.SeqIO.parse(handle, "fasta")])

def parse_genome(genome):
    scaffolds = []
    lengths = []
    for scaffold in genome.keys():
        scaffolds.append(scaffold)
        lengths.append(len(genome[scaffold]))
    return scaffolds, [length/float(sum(lengths)) for length in lengths], sum(lengths)

def weighted_choice(choices, weights):
    totals = []
    running_total = 0
    for w in weights:
        running_total += w
        totals.append(running_total)
    rnd = random.random() * running_total
    for choice, total in zip(choices, totals):
        if rnd < total:
            return choice

def gc_content(seq):
    seqlen = len(seq)
    gc = sum([1 if base in ['G', 'g', 'C', 'c'] else 0 for base in seq])
    return gc/float(seqlen)

def revcomp(seq):
    comp = {'a':'t', 'A':'T', 't':'a', 'T':'A', 'g':'c', 'G':'C', 'c':'g', 'C':'G'}
    return ''.join([comp[base] for base in seq[::-1]])

def amplify_frag(frag_seq, gc_to_bias):
    gc = find_nearest(compute_gcc(frag_seq), gc_to_bias.keys())
    return random.random() < gc_to_bias[gc]

def compute_gcc(seq):
    """ Return % GC content of string """
    gc = 0
    for base in seq:
        gc += 1 if base in ['G', 'g', 'C', 'c'] else 0
    return float(gc)/len(seq)

def find_nearest(value, list):
    """ Find nearest value of item in list """
    diffs = [abs(value - item) for item in list]
    index = diffs.index(min(diffs))
    return list[index]

def break_dna(dinuc, dinuc_to_bias):
    """
    Fragmentation bias
    Ref:
        http://www.nature.com/srep/2014/140331/srep04532/pdf/srep04532.pdf
        doi:10.1038/srep04532
    """
    return random.random() < dinuc_to_bias[str(dinuc).upper()]

def mutate_read(read, error_model, error_rate):
	""" Return read with error according to error model and error rate"""
	if error_model == 'illumina':
		return illumina_error(read)
	elif error_model == 'uniform':
		return uniform_error(read, error_rate)
	else:
		sys.exit('Unknown error model: %s' % error_model)

def illumina_error(read):
	"""
	Return read with Illumina sequencing error
	Reference: http://genomebiology.com/content/pdf/gb-2009-10-2-r23.pdf
	"""
	read_error = ''
	for i, base in enumerate(read):
		# error
		p = (3e-3 + 3.3e-8 * (i+1)**4.0)/100.00 
		if random.random() < p:
			if random.random() < 0.8: # 80% of errors are substitutions
				read_error += random.choice(['A','T','C','G'])
			elif random.random() < 0.5: # 10% of errors are insertions
				read_error += random.choice(['A','T','C','G'])
				read_error += base
			else: # 10% of errors are deletions
				continue
		else: # no error
			read_error += base
	return read_error

def uniform_error(read, error_rate):
	"""
	Return read with uniform sequencing error
	Reference: http://genomebiology.com/content/pdf/gb-2009-10-2-r23.pdf
	"""
	read_error = ''
	for i, base in enumerate(read):
		if random.random() < error_rate:
			if random.random() < 0.8: # 80% of errors are substitutions
				read_error += random.choice(['A','T','C','G'])
			elif random.random() < 0.5: # 10% of errors are insertions
				read_error += random.choice(['A','T','C','G'])
				read_error += base
			else: # 10% of errors are deletions
				continue
		else: # no error
			read_error += base
	return read_error












