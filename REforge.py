#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
from os.path import basename
import io
import sys
import logging
import time
import subprocess
import rpy2
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from Bio import Phylo, SeqIO
from Scoring_utilities import dollo_parsimony

__description__ = "Analyses transcription factor motifs with respect to their differences in binding within a phylogeny of CRM sequences and the associated phenotype"
__default_windowsize__ = 500
windowsize = __default_windowsize__
scorefile = "scores"

def __get_arguments():
	global scorefile, windowsize
	import argparse
	app = argparse.ArgumentParser(description=__description__)
	app.add_argument("treefile", 		type=str, help="phylogenetic tree given in newick format")
	app.add_argument("motiffile",		type=str, help="wtmx file containing all transcription factor motifs")
	app.add_argument("lossfile", 		type=str, help="file listing trait loss species; one line per species")
	app.add_argument("elementfile", 	type=str, help="file listing putative elements; one line per fasta file (.fa or .fasta)")

	app.add_argument("--scorefile", 		type=str, help="name of score file")
	app.add_argument("--windowsize", "-w",	type=str, help="windowsize to be used for element scoring")
	app.add_argument("--background", "-bg",	type=str, help="background file - passed into the score function")
	app.add_argument("--scrCrrIter", 		type=int, help="corrects stubb score with average score of <number> of shuffled sequences. Default is 10", default=10) 
	app.add_argument("--no_ancestral_filter", 	action="store_true", help="elements not are filtered if the ancestral element scores <= 0")
	app.add_argument("--no_branch_filter", 	action="store_true", help="branches are not filtered if start and end score <= filter_branch_threshold") 
	app.add_argument("--no_fixed_TP", 		action="store_true", help="Stubb's transition probablities are not fixed on the ancestral sequence") 

	app.add_argument("--filter_branch_threshold", 	type=float, help="Default is 0")
	# app.add_argument("--filter_branches", 		type=str, help="like --filter_branch_threshold but with motif specific branch thresholds")
	app.add_argument("--filter_GC_change",		type=float, help="if set, branches with a GC content change above the threshold are filtered")
	app.add_argument("--filter_length_change", 	type=float, help="if set, branches with a relative length change above the threshold are filtered")

	app.add_argument("--alignment", type=str)	

	app.add_argument("--verbose", "-v", action="store_true")
	app.add_argument("--debug", "-d", action="store_true")

	args = app.parse_args()
	if args.no_branch_filter:			assert args.filter_branch_threshold is None
	if args.windowsize is not None:		windowsize = args.windowsize
	if args.scorefile is not None: 		scorefile = args.scorefile
	return args

def create_branch_scoring_job_list(elements, tf, windowsize, args, filter_branch_threshold=None, background=None, stubb=False, append_to="jobFile"):
	with open(append_to, "a") as job_file:
		for cne in elements:
			arguments = ""
			if args.no_fixed_TP:						arguments += "--no_fixed_TP "
			if background is not None:					arguments += "-bg={0} ".format(background)
			# if filter_branches is not None:				arguments += "--filter_branches {0} ".format(filter_branches)
			if filter_branch_threshold is not None:		arguments += "--filter_branch_threshold {0} ".format(filter_branch_threshold)
			if args.filter_GC_change is not None:		arguments += "--filter_GC_change {0} ".format(args.filter_GC_change)
			if args.filter_length_change is not None:	arguments += "--filter_length_change {0} ".format(args.filter_length_change)
			if args.scrCrrIter is not None:				arguments += "--scrCrrIter {0} ".format(args.scrCrrIter)
			if args.no_ancestral_filter:				arguments += "--no_ancestral_filter "
			if args.no_branch_filter:					arguments += "--no_branch_filter "
			if args.alignment:							arguments += "--alignment {0} ".format(args.alignment)
			if args.debug:								arguments += "-d "
			if args.verbose:							arguments += "-v "
			job_file.write("echo -en {0}'\t'\"{1}\"'\t'; REforge_branch_scoring.py {0} {3} \"{1}\" {2} {4}; echo '<'\n".format(tf, cne, windowsize, args.treefile, arguments, ))

def __analyse_sequences():
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

	if os.path.exists(args.elementfile):
		elements = [line.strip() for line in open(args.elementfile, "r")]
	else:
		elements = args.elementfile.split(',')
	
	#########################
	# create a job file per transcription factor, for computing branch scores (eventually by extracting scores from the simulation) for every CNE and running them

	subprocess.check_call("> alljobs_%s"%scorefile, shell=True)													# (re-)create a new joblist file
	subprocess.check_call(["chmod", "+x", "alljobs_%s"%scorefile])
	allJobsFile = "alljobs_%s"%scorefile
	if not args.no_fixed_TP:
		subprocess.check_call(["mkdir", "-p", "tempDir"])

	if not os.path.exists(scorefile):	subprocess.check_call(["touch", scorefile])
	create_branch_scoring_job_list(elements, args.motiffile, windowsize, args, 
			filter_branch_threshold=args.filter_branch_threshold, #filter_branches=args.filter_branches,
			background=args.background, stubb=True, append_to=allJobsFile)

if __name__ == "__main__":	
	sys.exit(__analyse_sequences())