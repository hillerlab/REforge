#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import print_function
import os
import sys
import logging
import subprocess
from tempfile import mkstemp
from Bio import Phylo, SeqIO
from Scoring_utilities import score_sequence_with_stubb, computeGCContent, roundToList, read_sequences_and_prune_tree



__description__ = "Helper functions for the motif analysis pipeline. compute branch scores on a complete phylogeny."

sequence_score_cache = {"" : 0}

def __get_arguments():
	import argparse

	app = argparse.ArgumentParser(description=__description__)
	app.add_argument("motiffile",			type=str, help="motif as position frequency matrix")
	app.add_argument("treefile", 			type=str, help="phylogenetic tree. If elements might use a subtree ancestral nodes must be unnamed")
	app.add_argument("element", 			type=str, help="path to fasta file containing element's sequence alignment")							
	# species sequence data can be given in a BerkeleyDB file 'ALI'. In this case the parameter specifies the element name
	app.add_argument("scoring_window_length", 	type=str, help="scoring window length")
	app.add_argument("--background", "-bg",	type=str, help="background file - passed into the score function")
	app.add_argument("--no_ancestral_filter", 		action="store_true", help="elements are filtered if the ancestral element scores <= 0")
	app.add_argument("--no_fixed_TP", 		action="store_true", help="fix stubb's transition probablities on ancestral sequence") 
	app.add_argument("--filter_branch_threshold", "-dt",	type=float, help="if set, branches are filtered if start and end node score below the given threshold")
	app.add_argument("--filter_branches", "-db",	type=str, help="like --filter_branch_threshold but with motif specific branch thresholds - deprecated")
	app.add_argument("--filter_GC_change", "-fgc",	type=float, help="if set, branches with a GC content change above the threshold are filtered")
	app.add_argument("--filter_length_change", "-fl",	type=float, help="if set, branches with a relative length change above the threshold are filtered")
	app.add_argument("--scrCrrIter", 		type=int, help="corrects stubb score with average score of <number> of shuffled sequences. Default is 10", default=10) 

	app.add_argument("--verbose", "-v",		action="store_true")
	app.add_argument("--debug", "-d",		action="store_true")

	args = app.parse_args()
	return args


def read_threshold_collection(threshold_file_path, motif_name):		# returns {"stubb":stubb_dict)
	""" deprecated """
	stubb_dict = {}
	if not os.path.exists(threshold_file_path + motif_name + "_quantil.csv"):
		logging.fatal("threshold collection %s%s_quantil.csv not found"%(threshold_file_path, motif_name))
	with open(threshold_file_path + motif_name + "_quantil.csv", "r") as f:
		for line in f:
			parts = line.split("\t")
			stubb_dict[(float(parts[0]), motif_name)] = float(parts[1])
	return {"stubb" : stubb_dict}


def traverse_tree(parent, wtmx_file, scoring_window_length, sequences, background="", filter_threshold=None, 
					filter_threshold_dict=None, filter_GC_change=None, filter_length_change=None, 
					fixProbs=None, scrCrrMthd='stubb', scrCrrIter=10, anc0filter=False):
	""" traverses the phylegenetic tree, scores every sequence and computes branch scores """
	# anc0Filter: should be activated only for the root node
	# fixProbs: either a file, a bash string or True 
	global sequence_score_cache
	logging.info("checking {} with Stubb, scrCrrMthd: {}".format(parent.name, str(scrCrrMthd)) )
	sequence = str(sequences[parent.name].seq).replace("-", "")
	try:
		if sequence == '':
			return (0, 0, 0)
		if sequence.replace('N','') == '':
			return (float('nan'), 0, 0)
		gc = computeGCContent(sequence=sequence)/100
		if filter_threshold_dict is not None:
			gc_rounded = roundToList( [x[0] for x in filter_threshold_dict.keys()] , gc )
		else:
			gc_rounded = None

		if sequence in sequence_score_cache.keys():
			scores = [sequence_score_cache[sequence]]
		else:
			if fixProbs:
				if os.path.exists(str(fixProbs)):	
					scores = score_sequence_with_stubb(wtmx_file, scoring_window_length, sequence=sequence, 
							background=background, gc_content=gc, scrCrrMthd=scrCrrMthd, 
							scrCrrIter=scrCrrIter, fitprobs_file=fixProbs)
				else:
					fd, seq_tmpfile = mkstemp(prefix=parent.name + "-", dir="tempDir", suffix='.tmp')
					os.write(fd, bytes(">%s\n"%parent.name, 'utf-8') )
					os.write(fd, bytes( sequence, 'utf-8') )
					scores = score_sequence_with_stubb(wtmx_file, scoring_window_length, seq_file=seq_tmpfile, 
							background=background, gc_content=gc, scrCrrMthd=scrCrrMthd, 
							scrCrrIter=scrCrrIter, keepFiles=True)
					for suffix in ['', '.dict', '.fen', '.prof', '.parameters']:		# removes all created files except <seq_tmpfile>.fitprobs
						os.remove(seq_tmpfile + suffix)
					fixProbs = seq_tmpfile + ".fitprobs"
			else:
				scores = score_sequence_with_stubb(wtmx_file, scoring_window_length, sequence=sequence, background=background,
						gc_content=gc, scrCrrMthd=scrCrrMthd, scrCrrIter=scrCrrIter)
			logging.debug(scores)

			sequence_score_cache[sequence] = scores[0]
		if len(scores) > 0:	parent_score = sum(scores) / len(scores)
		else:  				parent_score = 0  											# in this case no sequence was given in the fasta file

		# if root does not contain binding sites and fixed transition probabilities are activated tree traversal can be skipped
		if anc0filter and parent_score <= 0: 
			return (parent_score, gc, len(sequence))

		for child in parent:
			child_score, child_gc, child_length = traverse_tree(child, wtmx_file, scoring_window_length, sequences,
						background, filter_threshold, filter_threshold_dict, filter_GC_change, filter_length_change, 
						fixProbs, scrCrrMthd=scrCrrMthd, scrCrrIter=scrCrrIter, anc0filter=False)

			motif_name = os.path.splitext(os.path.basename(wtmx_file))[0]
			# branch score filter step
			if not any([filter_threshold_dict is not None and \
						parent_score < filter_threshold_dict[(gc_rounded, motif_name)] and \
						child_score < filter_threshold_dict[(gc_rounded, motif_name)],
					filter_threshold is not None and parent_score < filter_threshold and child_score < filter_threshold,
					filter_GC_change is not None and abs(child_gc - gc) > filter_GC_change,
					filter_length_change is not None and abs(child_length - length) < filter_length_change*(child_length+length)/2 ]): 
				print("%s>%s:%f\t"%(parent, child, child_score  - parent_score), end="")

			logging.debug("%s>%s:%f\t(child_score: %f\tparent_score: %f)"%(parent, child, child_score - parent_score, child_score, parent_score))
	except:	raise
	finally:
		try:
			os.remove( seq_tmpfile + ".fitprobs")
		except (NameError, FileNotFoundError):	pass

	return (parent_score, gc, len(sequence))


def __score_branches():
	"""the main function"""
	import io

	# get the arguments and handle --verbose logging.
	args = __get_arguments()
	if args.debug:		logging.basicConfig(level=logging.DEBUG)
	elif args.verbose:	logging.basicConfig(level=logging.INFO)
	else:				logging.basicConfig(level=logging.WARNING)

	logging.info("branch scoring of {0} with {1}".format(args.element, args.motiffile))
	

	# read the input files
	phylo_tree, sequences = read_sequences_and_prune_tree(args.element, args.treefile)

	if args.filter_branches:
		filter_threshold_dict = read_threshold_collection(args.filter_branches, os.path.splitext(os.path.basename(args.motiffile))[0])["stubb"]
	else:
		filter_threshold_dict = None

	logging.info("Branch scoring: use %s as background and %s as window size"%(args.background, args.scoring_window_length) )
	scrCrrMthd = "stubb" if args.scrCrrIter > 0 else None
	traverse_tree(phylo_tree.root, args.motiffile, float(args.scoring_window_length), sequences, args.background, 
					filter_threshold=args.filter_branch_threshold, filter_threshold_dict=filter_threshold_dict,
					filter_GC_change=args.filter_GC_change, filter_length_change=args.filter_length_change,
					fixProbs=not args.no_fixed_TP, scrCrrMthd=scrCrrMthd, scrCrrIter=args.scrCrrIter,
					anc0filter=not args.no_ancestral_filter)

	logging.info(computeGCContent.cache_info())
	return 0


if __name__ == "__main__":
	sys.exit(__score_branches())
