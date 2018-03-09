#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
import sys
import logging
import math
import subprocess
import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from Bio import SeqIO 
from REforge_util import read_sequences_and_prune_tree, dollo_parsimony

__description__ = "Analyses transcription factor motifs with respect to their differences in binding within a phylogeny of CRM sequences and the associated phenotype"
suffix = ""									# suffix ending for newly generated files

def __get_arguments():
	global suffix
	import argparse
	app = argparse.ArgumentParser(description=__description__)

	app.add_argument("treefile", 		type=str, help="phylogenetic tree given in newick format")
	app.add_argument("motiffile",		type=str, help="wtmx file containing all transcription factor motifs")
	app.add_argument("lossfile", 		type=str, help="file either in special loss format (.ls) or simple list")
	app.add_argument("elementfile", 	type=str, help="file with names of putative elements; one line per element")

	app.add_argument("--add_suffix", 	type=str, help="add this suffix to files")
	app.add_argument("--filterspecies", type=str, help="skip species")
	app.add_argument("--elements",		type=str, help="analyze only the elements specified in this file")

	app.add_argument("--verbose", "-v", action="store_true")
	app.add_argument("--debug", "-d", action="store_true")

	args = app.parse_args()
	if args.add_suffix is not None: 	suffix = args.add_suffix
	return args

def test_significance(s, p):
	""" computes significance of a positive Pearson correlation between s and (binary) p """
	# cohensd		= robjects.r('function (x, y) { cd <- (mean(x) - mean(y)) / sd(y) }')
	# maxZscore	= robjects.r('function (x, y) { cd <- (max(x) - mean(y)) / sd(y) }')

	r_s = robjects.FloatVector(s)
	r_p = robjects.FloatVector(p)
	data = robjects.DataFrame({"Val":r_s, "P":r_p})					# is ordered alphabetically, ie data[0] corresponds to P
	data0 = robjects.r['subset'](data, data.rx('P').ro == 0)
	data1 = robjects.r['subset'](data, data.rx('P').ro == 1)
	# t 	= robjects.r['t.test'](robjects.r('Val ~ P'), alternative="less", data=data )
	Pearson		= robjects.r['cor.test'](data.rx2('P'), data.rx2('Val'), method="pearson", exact=True, alternative="greater")
	# Spearman	= robjects.r['cor.test'](data.rx2('P'), data.rx2('Val'), method="spearman", exact=True, alternative="greater")
	# try:
	# 	w 	= robjects.r['wilcox.test'](robjects.r('Val ~ P'), data=data, alternative="less", **{'conf.int': True})
	# except rpy2.rinterface.RRuntimeError as err:
	# 	logging.warning(err)
	# 	w 	= None

	# if len(data0.rx2('Val')) > len(data1.rx2('Val')):
	# 	logging.fatal("less background than foreground scores - background scores are used as variance estimation for Cohen's D")
	# cd 	= cohensd(data0.rx2('Val'), data1.rx2('Val'))
	# zs  = maxZscore(data0.rx2('Val'), data1.rx2('Val'))
	# ave = (Pearson.rx2('estimate')[0] + Spearman.rx2('estimate')[0]) / 2
	# df = len(p) - 2
	# statistic = robjects.r['c']( t = math.sqrt(df) * ave/math.sqrt(1 - ave**2) )
	# PearsonSpearmanP = 1 - robjects.r['pt'](statistic, df)[0]

	# texts = []
	# for p in [PearsonSpearmanP, Pearson.rx2('p.value')[0], Spearman.rx2('p.value')[0], w.rx2('p.value')[0] if w is not None else 'NA', t.rx2('p.value')[0]]:
	# 	if p == 'NA':	texts.append("NA")
	# 	else:	 		texts.append("%e"%p)
	# if w is None:	texts.append("NA")
	# else:			texts.append("%e"%w.rx2('conf.int').rx(2)[0])
	# texts.extend([ "%e"%t.rx2('conf.int').rx(2)[0], "%f"%cd[0], "%f"%zs[0] ])
	ret_val = str(Pearson.rx2('p.value')[0])
	# ret_val = "\t".join(texts)
	return ret_val
	
def compute_and_output_significance(args, scoreData, elements, more_data=None):
	""" creates computes significance for every element and creates summarizing file """
	full_suffix = suffix
	if args.elements: full_suffix += '_' + os.path.basename(args.elements)

	with open("significant_elements%s"%full_suffix, "w") as summary:
		summary.write("element\tPearson\t#pos\t#neg\n")
		for cne in elements:
			if cne+"+" in scoreData and cne+"-" in scoreData:
				phenotypes = [1]*len(scoreData[cne+"+"]) + [0]*len(scoreData[cne+"-"])
				scores = scoreData[cne+"+"] + scoreData[cne+"-"]
				if phenotypes.count(0) > 1 and phenotypes.count(1) > 1:
					try:
						pvalue = test_significance(scores, phenotypes)
					except Exception as err:
						logging.warning("Error in test_significance of %s"%cne)
						logging.warning("scores: %s"%(str(scores)))
						logging.warning("phenotypes: %s"%(str(phenotypes)))
						raise err
					summary.write("{0:15}\t{1}\t{2}\t{3}\n".format(cne, pvalue, phenotypes.count(0), phenotypes.count(1)))
				else:
					summary.write("{0:15}\tNA\t{2}\t{3}\n".format(cne, phenotypes.count(0), phenotypes.count(1)) )
			else:
				logging.info("scoreData for transcription factor %s is missing"%cne)


def __analyse_sequences():
	import functools

	args = __get_arguments()
	if args.debug:
		logging.basicConfig(level=logging.DEBUG)
		logging.debug("logging level is DEBUG")
	elif args.verbose:
		logging.basicConfig(level=logging.INFO)	
		logging.info("logging level is INFO")
	else:
		logging.basicConfig(level=logging.WARNING)
		logging.debug("logging level is WARNING")

	
	elements = [line.strip() for line in open(args.elementfile, "r")]
	transcription_factors = [args.motiffile]

	#########################
	# read scores
	logging.info("reading scorefile")
	scoreFileContent = {}
	if args.elements:
		logging.info("cropping scorefile down to element list")
		sf = subprocess.check_output("grep -f %s scores%s"%(args.elements, suffix), shell=True).decode('utf-8').strip().split('\n')
	else:
		sf = open("scores%s"%suffix)

	counter=0
	c_nan = 0
	for line in sf:
		counter += 1
		parts = line.split()
		if parts[-1] != '<':	logging.fatal("Malformed scores%s file - line does not end with '<': %s"%(suffix, line))
		tf, cne = parts[0], parts[1]
		if (tf, cne) in scoreFileContent:	logging.fatal("Malformed scores%s file - two entries for TF-CNE pair: %s"%(suffix, line))	
		if cne in elements:
			scoreFileContent[tf,cne] = {}
			try:
				for part in parts[2:-1]:
					branch, score = part.split(':')
					if math.isfinite(float(score)):
						scoreFileContent[(tf,cne)][branch] = float(score)
					else:
						c_nan += 1
			except Exception as err:
				logging.fatal("Error occured while parsing scores%s at line '%s'"%(suffix, line))
				raise err
	if not args.elements: sf.close()
	logging.info("finished reading scorefile")
	if c_nan > 0:	logging.warning("%d NaN or Inf values have been skipped"%c_nan)
	logging.info("%d lines read"%counter)

	#########################
	# retrieve trait loss branches
	with open(args.lossfile) as lf:
		losses = [line.strip() for line in lf]
	dismiss_species = args.filterspecies.split(",") if args.filterspecies else []
	sub_phylos = {}
	for cne in elements:
		pruned_tree, _ = read_sequences_and_prune_tree(cne, args.treefile, dismiss_species)
		sub_phylos[cne] = dollo_parsimony(pruned_tree, losses)			# annotate with trait
	logging.info("trait loss information read and computed")

	#########################
	# read, classify and analyse scores
	scoreData	= {key+"-": [] for key in elements}
	scoreData.update({key+"+": [] for key in elements})

	# skip_cne	= True
	cache = {}
	ch, cm = 0,0 
	
	for cne in elements:
		branch_ends = [n.name for n in sub_phylos[cne].get_nonterminals()] + [n.name for n in sub_phylos[cne].get_terminals()]
		branch_ends.remove(sub_phylos[cne].root.name)
		for branch_end in branch_ends:
			branch = ([sub_phylos[cne].root] + sub_phylos[cne].get_path(branch_end))[-2].name +'>' + branch_end
			if (str(sub_phylos[cne]), branch_end) in cache:
				trait = cache[(str(sub_phylos[cne]), branch_end)]
				ch += 1
			else:
				trait = sub_phylos[cne].find_any(branch_end).trait
				cache[(str(sub_phylos[cne]), branch_end)] = trait
				cm += 1
			for tf in transcription_factors:
				if branch in scoreFileContent[(tf,cne)]:
					score = scoreFileContent[(tf,cne)][branch]
					if trait:	scoreData[cne + "+"].append(score)
					else: 		scoreData[cne + "-"].append(score)

	logging.debug("Hits: %d; Misses: %d"%(ch, cm))

	#########################
	# assess phenotype-divergence association
	compute_and_output_significance(args, scoreData, elements)

if __name__ == "__main__":	
	sys.exit(__analyse_sequences())