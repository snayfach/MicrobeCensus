# MicrobeCensus - estimation of average genome size from shotgun sequence data
# Copyright (C) 2013-2014 Stephen Nayfach
# Freely distributed under the GNU General Public License (GPLv3)

"""
Usage:
	python optimize_parameters.py

Description: 
	find the "optimal", read-length specific mapping parameters for each marker gene
	optimal parameters are those which result in the most accurate estimates of AGS

"""

# Libraries
import sys, os, numpy
from training_functions import *

# Arguments
threads = int(sys.argv[1])
xfolds  = int(sys.argv[2])

# Filepaths to input/output data
main_dir  = os.path.dirname(os.path.dirname(os.path.realpath(sys.argv[0])))
src_dir   = os.path.join(main_dir, 'src')
train_dir = os.path.join(main_dir, 'training')
input_dir = os.path.join(train_dir, 'input')
inter_dir = os.path.join(train_dir, 'intermediate')
out_dir   = os.path.join(train_dir, 'output')

genomes_dir = os.path.join(input_dir, 'genomes')
reads_dir   = os.path.join(inter_dir, 'reads')
hits_dir    = os.path.join(inter_dir, 'hits')

pars_path   = os.path.join(out_dir, 'pars.map')
coeffs_path = os.path.join(out_dir, 'coefficients.map')
preds_path  = os.path.join(out_dir, 'training_preds.map')

# Check if input directories exist
for dir in [src_dir, train_dir, input_dir, inter_dir, out_dir, genomes_dir, reads_dir, hits_dir]:
	if not os.path.isdir(dir):
		os.mkdir(dir)

# Store genome & library sizes
print "Storing genome sizes and library sizes..."
genome2size = genome_sizes(genomes_dir)
library2size = library_sizes(reads_dir)

# Store classification rates
print "Storing classification rates..."
class_rates = store_rates(hits_dir, library2size)

# Build argument list
args_list = []
for read_length in class_rates:
	for fam in class_rates[read_length]:
		for min_score, max_pid, aln_cov in class_rates[read_length][fam]:
			for rate_type in ['rate_hits', 'rate_aln', 'rate_cov']:
				arguments = {}
				arguments['pars'] = (read_length, fam, min_score, max_pid, aln_cov, rate_type)
				arguments['xfolds'] = xfolds
				arguments['genome_names'] = class_rates[read_length][fam][(min_score, max_pid, aln_cov)]['genome_names']
				arguments['rates'] = class_rates[read_length][fam][(min_score, max_pid, aln_cov)][rate_type]
				arguments['genome2size'] = genome2size
				args_list.append(arguments)

# Compute x-validation error for each set of parameters
# and find optimal parameters

print "Evaluating performance for %s parameters..." % len(args_list)
xval_error = parallel_return_function(xvalidation, args_list, threads)

print "Finding optimal parameters..."
opt_pars = find_opt_pars(xval_error)
#	write results
f_out = open(pars_path, 'w')
header = ['gene_fam', 'read_length', 'aln_cov', 'max_pid', 'min_score', 'aln_stat']
f_out.write('\t'.join(header)+'\n')
for read_length, fam in sorted(opt_pars.keys()):
	min_score, max_pid, aln_cov, rate_type = opt_pars[(read_length, fam)]['pars']
	aln_stat = 'hits' if rate_type == 'rate_hits' else 'aln' if rate_type == 'rate_aln' else 'cov'
	record = [ str(x) for x in [fam, read_length, aln_cov, max_pid, min_score, aln_stat] ]
	f_out.write('\t'.join(record)+'\n')

print "Finding proportionality constants..."
prop_consts = {}
for read_length, fam in sorted(opt_pars.keys()):
	min_score, max_pid, aln_cov, rate_type = opt_pars[(read_length, fam)]['pars']
	genome_names = class_rates[read_length][fam][(min_score, max_pid, aln_cov)]['genome_names']
	rates = class_rates[read_length][fam][min_score, max_pid, aln_cov][rate_type]
	prop_const = estimate_proportionality_constant(genome_names, rates, genome2size)
	prop_consts[read_length, fam] = prop_const
#	write results
f_out = open(coeffs_path, 'w')
for read_length, fam in sorted(prop_consts.keys()):
	prop_const = prop_consts[(read_length, fam)]
	record = [read_length + '_' + fam, str(prop_const)]
	f_out.write('\t'.join(record)+'\n')

print "Getting per-gene AGS predictions..."
training_ests = {}
for read_length, fam in sorted(opt_pars.keys()):
	min_score, max_pid, aln_cov, rate_type = opt_pars[(read_length, fam)]['pars']
	prop_const = prop_consts[(read_length, fam)]
	genome_names = class_rates[read_length][fam][(min_score, max_pid, aln_cov)]['genome_names']
	rates = class_rates[read_length][fam][min_score, max_pid, aln_cov][rate_type]
	for genome_name, rate in zip(genome_names, rates):
		true_ags = genome2size[genome_name]
		est_ags = prop_const/rate if rate > 0 else 'NA'
		training_ests[(read_length, fam, genome_name)] = true_ags, est_ags
#	write results
f_out = open(preds_path, 'w')
header = ['read_length', 'fam', 'genome_name', 'true_ags', 'est_ags']
f_out.write('\t'.join(header)+'\n')
for read_length, fam, genome_name in training_ests:
	true_ags, est_ags = training_ests[(read_length, fam, genome_name)]
	record = [ str(x) for x in [read_length, fam, genome_name, true_ags, est_ags] ]
	f_out.write('\t'.join(record)+'\n')
