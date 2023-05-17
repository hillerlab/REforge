#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import logging
import numpy as np
import os
import random
import re
import subprocess
import sys
from tempfile import mkstemp
from functools import lru_cache
from os import path
from Bio import Phylo, SeqIO, SeqUtils

__description__ = "helper for the REforge and TFforge pipeline to score elements"

def dollo_parsimony(phylo_tree, traitLossSpecies):
	""" simple dollo parsimony implementation
	Given a set of species that don't have a trait (in our hash), we do a bottom up pass through the tree.
	If two sister species have lost the trait, then the ancestor of both also lost it. Otherwise, the ancestor has the trait. """
	for node in phylo_tree.get_terminals():
		if node.name in traitLossSpecies:
			node.trait = False
		else:
			node.trait = True

	for node in phylo_tree.get_nonterminals(order='postorder'):
		child_traits = [child.trait for child in node]
		node.trait = bool(sum(child_traits))
	return phylo_tree


def read_sequences_and_prune_tree(fastafile, treefile, dismiss_species=[], alignment=None):
	""" reads the fastafile and prunes the tree to all species that occur  """
	import io
	sequences = {}
	try:
		if alignment is None:
			if not (fastafile.endswith('.fa') or fastafile.endswith('.fasta')): logging.warning("%s specified as element seems to be not in fasta-format."%fastafile)
			sequences = SeqIO.to_dict(SeqIO.parse(fastafile, "fasta"))
		else:
			sequencedata = subprocess.check_output("ReadBDB.perl %s '%s' -brief"%(alignment, fastafile), shell=True)
			sequences = SeqIO.to_dict(SeqIO.parse(io.StringIO(sequencedata.decode("utf-8")), "fasta"))
	except Exception as err:
		logging.fatal(err)
		
	# prune and read tree
	if sequences == {}:
		logging.fatal("no sequences given for tree pruning - empty tree")
		phylo_tree = Phylo.read(io.StringIO("()"), "newick")
	else:
		logging.debug("tree_doctor %s -antN -P %s"%(treefile, ",".join(sequences.keys())))
		treedata = subprocess.check_output("tree_doctor %s -antN -P %s"%(treefile, ",".join(sequences.keys() - set(dismiss_species))), shell=True)
		phylo_tree = Phylo.read(io.StringIO(treedata.decode("utf-8")), "newick")
	return phylo_tree, sequences


@lru_cache(maxsize=1000)
def computeGCContent(seq_file=None, sequence=None):
	""" computes the GC-content of a given sequence or sequence file. Returns the GC-content in % """
	if (seq_file is None) == (sequence is None):
		raise Exception("Error in computeGCContent: Either seq_file or sequence must be specified")
	if seq_file is not None:
		gc_contents = [SeqUtils.GC123(s.seq)[0] for s in SeqIO.parse(seq_file, format="fasta")]				# use GC123 instead of GC to cope with dashes
		if len(gc_contents) > 1:
			logging.debug("gc_content is averaged over all sequences found in %s"%seq_file)
		return sum(gc_contents)/len(gc_contents)
	else:
		return SeqUtils.GC123(sequence)[0]

def generate_random_sequence(length, prob_dist=[0.25, 0.25, 0.25, 0.25]):
	""" generates a DNA sequence given length  according to the given nucleotide distribution """
	seq = ''
	for p in range(length):
		n = random.random()
		if n < prob_dist[0]: seq += 'A'
		elif n < sum(prob_dist[0:2]): seq += 'C'
		elif n < sum(prob_dist[0:3]): seq += 'G'
		else: seq += 'T' 
	return seq

def shuffle_sequence(seq):
	""" shuffles a sequence given as string """
	l = list(seq)
	random.shuffle(l)
	return ''.join(l)

def roundToList(elements, item, returnIndex=False):
	""" rounds item to the its closest number that the list elements contain """
	np_el = np.array(sorted(elements))
	if returnIndex:	return (np.abs(np_el - item)).argmin()
	else:			return np_el[ (np.abs(np_el - item)).argmin() ]

def find_best_background(folder, gc_content=None, seq_file=None, sequence=None):
	""" searches a folder the subfolder of format 'gc0.xx' which matches best the given GC content or the GC content of the given sequence.
		The folder is expected to contain an equally named fasta file. """
	if gc_content is None and seq_file is None and sequence is None:
		raise Exception("Error in find_best_background: One of gc_content, seq_file or sequence must be specified")
	pattern = re.compile("(?<=gc)[0-9].[0-9]+")
	if gc_content is None:	gc = computeGCContent(seq_file=seq_file, sequence=sequence) / 100
	else:					gc = gc_content
	if not os.path.isabs(folder):	folder = os.getcwd() + '/' + folder
	if not os.path.exists(folder):
		logging.fatal("background %s not found"%folder)
	backgrounds = sorted( next( os.walk(folder) )[1] )
	bkgd_gcs = [re.search(pattern, subfolder) for subfolder in backgrounds]
	if bkgd_gcs == [None] * len(bkgd_gcs):
		logging.fatal("No folder name of the pattern 'gc0.xx or gc1.00' was found in the background folder %s"%folder)
	
	bkgd_ind = roundToList([ float(x.group(0)) for x in bkgd_gcs], gc, returnIndex=True)
	logging.debug("background {0}/{1}/{1}.fa will be used".format(folder, backgrounds[bkgd_ind]))
	return "{0}/{1}/{1}.fa".format(folder, backgrounds[bkgd_ind])

@lru_cache(maxsize=1000)
def score_sequence_with_stubb(wtmxFile, window, seq_file=None, sequence=None, background=None, gc_content=None, scrCrrMthd=None, scrCrrIter=10, fitprobs_file=None, keepFiles=False, arguments=""):
	if (seq_file is None) == (sequence is None):		raise Exception("Error in score_sequence_with_stubb: Either seq_file or sequence must be specified")
	if sequence is None and not path.exists(seq_file):	raise Exception("Error in score_sequence_with_stubb: sequence file %s not found"%seq_file)
	if glob.glob(wtmxFile) == []:						raise Exception("Error in score_sequence_with_stubb: wtmx file %s not found"%wtmxFile)
	if keepFiles and not seq_file:						raise Exception("Error in score_sequence_with_stubb: Files can only be kept with sequence given as file")
	assert scrCrrMthd in [None, 'stubb']
	
	# get background file
	if background is not None:
		if background.endswith(".fa"):	bkgd_file = background
		else:							bkgd_file = find_best_background(background, gc_content, seq_file=seq_file, sequence=sequence)
		arguments += ' -b ' + bkgd_file
	else:
		bkgd_file = None
		
		
	if scrCrrMthd is not None:
		# compute score correction value via application of Stubb on shuffled sequences
		logging.debug("Iterative score correction")
		if sequence is None:
			file_sequences = [str(s.seq) for s in SeqIO.parse(seq_file, format="fasta")]
			sequence_template = ''.join(file_sequences).replace('-','')
		else:
			sequence_template = sequence.replace('-','')
		random.seed(sequence_template)
		sequence_template = sequence_template.replace('N','')
		gc = computeGCContent(seq_file=seq_file, sequence=sequence) / 100 if gc_content is None else gc_content
		
		r_scores = []
		for x in range(scrCrrIter):
			r_seq	= shuffle_sequence(sequence_template)
			r_score = score_sequence_with_stubb(wtmxFile, window, sequence=r_seq, background=bkgd_file, fitprobs_file=fitprobs_file, scrCrrMthd=None)
			assert len(r_score) == 1 																# number of scores correspond to number of sequences
			r_scores.append( r_score[0] )
		r_score = sum(r_scores)/len(r_scores)
		
	# compute Stubb score
	if sequence is not None:
		if fitprobs_file:	call = "stubb_fixedprobs_noPseudoCount %s %s %s 1 %s %s -brief -seq"%(sequence, wtmxFile, window, fitprobs_file, arguments)
		else:				call = "stubb_noPseudoCount %s %s %s 1 %s -brief -seq"%(sequence, wtmxFile, window, arguments)
		logging.debug(call)
		try:
			output = subprocess.check_output(call, shell=True, executable='/bin/bash')
			output = list(map(float, output.strip().split()))
			scores = [output[1]]
		except Exception as err:
			logging.warning("Error in score_sequence_with_stubb: %s"%err)
			# Segmentation Fault occurs if number of N in window > length/2 where length is min(sequence length, window size)
			return [float('nan')]
	else:
		if fitprobs_file:	call = "stubb_fixedprobs_noPseudoCount %s %s %s 1 %s %s"%(seq_file, wtmxFile, window, fitprobs_file, arguments)
		else:				call = "stubb_noPseudoCount %s %s %s 1 %s"%(seq_file, wtmxFile, window, arguments)
		if keepFiles:
			logging.debug(call)
			try:
				subprocess.check_call(call, shell=True, executable='/bin/bash')
			except subprocess.CalledProcessError as err:
				logging.warning("Error occured while scoring sequence file %s :\n %s"%(seq_file, err))
				return [float('nan')]
			try:
				with open(seq_file + ".fen") as f:
					scores = []
					for x in f:
						lps = x.split()
						if lps[0] == '0':	scores.append(float(lps[1]))
						else:				scores[-1] = max(scores[-1], float(lps[1]))
					if scores == []:	
						if scrCrrMthd: assert r_score == 0
						return [0] 									# sequence was too short to be scored
			except ValueError:
				raise Exception("Error in %s"%(call))
		else:
			call += ' -brief'
			logging.debug(call)
			output = subprocess.check_output(call, shell=True, executable='/bin/bash')
			scores = [list(map(float, output.strip().split()))[1]]
			
	if scrCrrMthd is not None:	return [s - r_score for s in scores]
	else:						return scores
