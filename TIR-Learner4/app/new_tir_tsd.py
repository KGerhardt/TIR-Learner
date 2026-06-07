import sys
import os

#from new_seq_reader import json_loader, bed_worker

import multiprocessing
import re

from pywfa import WavefrontAligner

import numpy as np

#from pywfa import cigartuples_to_str

class tsd_tir_checker:
	def __init__(self):
		#Because CNN is so expensive, it makes sense to check all TIRvish + GRF results for validity
		self.invert = {}
	
		self.rc_trans = str.maketrans('ATCG', 'TAGC')
	
		self.tir_pattern_dict = None
		self.no_tir_motif = set(["DTX", "NonTIR"])
		dta_pattern = re.compile(r'[CT]A[AG][ATGC]G')
		dtc_pattern = re.compile(r'CACT[AG]')
		dth_pattern = re.compile(r'G[GA][GC]C')
		dtm_pattern = re.compile(r'[GC]')
		dtt_pattern = re.compile(r'CT[ATCG][ATCG]CTC[ATCG][ATCG]T')
		dte_pattern = re.compile(r'GG[ATCG][AG][AC]')
		dtr_pattern = re.compile(r'CAC[AT]ATG')
		dtp_pattern = re.compile(r'CA[ATGC][AG]G')

		self.tir_pattern_dict = {
						'DTA':dta_pattern,
						'DTC':dtc_pattern,
						'DTH':dth_pattern,
						'DTM':dtm_pattern,
						'DTT':dtt_pattern,
						'DTE':dte_pattern,
						'DTR':dtr_pattern,
						'DTP':dtp_pattern}
		
		self.tsd_min_size_by_family = {
										"DTA": 8, #[8],
										"DTC": 2, #[3, 2],
										"DTH": 3, #[3],
										"DTM": 7, #[10, 9, 8, 7],
										"DTT": 2,
										"DTE": 2,
										"DTR": 2,
										"DTP": 2,
										"DTX": 2} #[2]}
		
		self.no_tsd_motif = set(['DTX', 'NonTIR', 'DTC', 'DTM','DTE', 'DTR', 'DTP'])
		
		self.tsd_pattern_dict = {'DTH': re.compile(r'T[TA]A'),
								'DTA': re.compile(r'TA')}
		
	def revcomp(self, seq):
		rc = seq[::-1].translate(self.rc_trans)
		return rc
		
	#This is used as a checker for the mere existence of a repeat in the left/right sequences above a similarity cutoff;
	#It is one half of the checks that were done in the original check_TIR_TSD, and the separation is deliberate
	#The point being, this check can be applied to TIRvish + GRF sequences BEFORE they are passed to more expensive 
	#steps of BLAST or CNN, then only the post-checks need applied afterwards
	
	#Functionally this is just for checking TIRs
	def wfa_align(self, left, right, min_size = 10, min_similarity = 0.8):
		aln = WavefrontAligner(left)
		result = aln(right)
		
		calc_array = np.array(result.cigartuples, dtype = np.int32)
		repeat_length_arr = np.cumsum(calc_array[:, 1])
		repeat_match_count = np.zeros(calc_array.shape[0], dtype = np.int32)
		match_indices = calc_array[:, 0] == 0
		repeat_match_count[match_indices] = calc_array[match_indices, 1]
		repeat_match_count = np.cumsum(repeat_match_count)
		
		#These are the similarity indices
		rolling_percent_sim = repeat_match_count / repeat_length_arr

		#We want the max of this that still ends with a match
		acceptable_percent_and_length = np.logical_and(rolling_percent_sim > min_similarity, repeat_length_arr >= min_size)
		
		has_repeat = False
		left_repeat_size = None
		right_repeat_size = None
		ref_starts_at = None
		query_starts_at = None
		tir_percent = None
		
		if np.any(acceptable_percent_and_length):
			has_repeat = True
			final_ok = np.where(acceptable_percent_and_length)[0][-1]
			while final_ok > 0:
				#This checks if the array ends with one or several mismatches of any kind and backtracks to the last match
				if acceptable_percent_and_length[final_ok - 1] > acceptable_percent_and_length[final_ok]:
					final_ok -= 1
				else:
					break
			
			#Total mismatches
			shared_mismatches = np.sum(calc_array[:final_ok, 1][calc_array[:final_ok, 0] == 8])
			
			#Matches + mismatches + insertions into ref / deletions from query
			left_repeat_size = repeat_match_count[final_ok] + shared_mismatches + np.sum(calc_array[:final_ok, 1][calc_array[:final_ok, 0] == 2])
			#matches + mismatches + insertions into query / deletions from ref
			right_repeat_size = repeat_match_count[final_ok] + shared_mismatches + np.sum(calc_array[:final_ok, 1][calc_array[:final_ok, 0] == 1])
			
			left_repeat_size  = int(left_repeat_size)
			right_repeat_size = int(right_repeat_size)
			
			as_dict = result.__dict__
			
			ref_starts_at = as_dict['pattern_start']
			query_starts_at = as_dict['text_start']
			
			tir_percent = round(100 * rolling_percent_sim[final_ok], 1)
		
		return has_repeat, left_repeat_size, right_repeat_size, ref_starts_at, query_starts_at, tir_percent
		
	def check_tir_conservation(self, tir_type, fwd, rev):
		is_conserved = False
		if tir_type in self.no_tir_motif:
			is_conserved = True
		else:
			matf = re.match(self.tir_pattern_dict[tir_type], fwd)
			matr = re.match(self.tir_pattern_dict[tir_type], rev)
			
			if matf or matr:
				is_conserved = True
				
		return is_conserved

	#This code needs to run the alignment, then proceed from the end of the left and start of the right to check
	#for the TSD starting at the end of the left seq and start of the right seq, proceed in opposite directions along both
	#The code also checks DTH and DTA TSDs for acceptable starts
	def check_tsd(self, left, right, tir_type, min_similarity = 0.8):
		if tir_type in self.tsd_min_size_by_family:
			min_ok_size = self.tsd_min_size_by_family[tir_type]
		else:
			min_ok_size = 2
			
		needs_conservation_check = tir_type in self.tsd_pattern_dict
		
		#max_mismatch = int((min_ok_size * (1-min_similarity)) - 0.5)
		
		has_tsd = False
		#maxl = min([len(left), len(right)])
		
		aln = WavefrontAligner(left, text_end_free=1, pattern_begin_free=1)
		
		winning_tsd = -1
		left_chars = None
		right_chars = None
		tsd_percent = None
		
		for i in range(min_ok_size, len(right) + 1):
			result = aln(right[:i])
			ciglen = len(result.cigartuples)
			#This is a full-length perfect match; couldn't be better
			if ciglen == 1:
				if result.cigartuples[0][0] == 0:
					if needs_conservation_check:
						if re.match(self.tsd_pattern_dict[tir_type], left[-i:]) or re.match(self.tsd_pattern_dict[tir_type], right[i:]):
							winning_tsd = i
							left_chars = i
							right_chars = i
							tsd_percent = 100.0
					else:
						winning_tsd = i
						left_chars = i
						right_chars = i
						tsd_percent = 100.0
			
			#this is a perfect match at a greater length than any yet checked in exactly the correct location; this is always a winner
			if ciglen == 2:
				if result.cigartuples[0][0] == 2 and result.cigartuples[1][0] == 0:
					if needs_conservation_check:
						if re.match(self.tsd_pattern_dict[tir_type], left[-i:]) or re.match(self.tsd_pattern_dict[tir_type], right[i:]):
							winning_tsd = i
							left_chars = i
							right_chars = i
							tsd_percent = 100.0
					else:
						winning_tsd = i
						left_chars = i
						right_chars = i
						tsd_percent = 100.0
			
			if ciglen > 2:
				#Only consider reassigning winning TSD if there's at least two extra hits gained
				if i > winning_tsd + 1:
					#Must be gappy against ref to first aligned char
					if result.cigartuples[0][0] == 2:
						match_sum = 0
						totl_size = 0
						tmp_left_chars = 0
						tmp_right_chars = 0
						for op, ct in result.cigartuples[1:]:
							#Regardless of match, mismatch, or gap, we consider it a an extension of length
							totl_size += ct
							#This is a match, consumes chars for both and adds to match sum
							if op == 0:
								match_sum += ct
								tmp_left_chars += ct
								tmp_right_chars += ct
							#This is a deletion, so it only consumes right chars
							if op == 1:
								tmp_right_chars += ct
							#This is an insertion, only consumes left chars
							if op == 2:
								tmp_left_chars += ct
							#This is a mismatch, consumes chars for both
							if op == 8:
								tmp_left_chars += ct
								tmp_right_chars += ct
						
						#the match contains a mismatch or gap, but is still good enough
						if match_sum > (totl_size * min_similarity):
							if needs_conservation_check:
								if re.match(self.tsd_pattern_dict[tir_type], left[-i:]) or re.match(self.tsd_pattern_dict[tir_type], right[i:]):
									winning_tsd = i
									left_chars = tmp_left_chars
									right_chars = tmp_right_chars
									tsd_percent = round(100 * match_sum / totl_size, 1)
							else:
								winning_tsd = i
								left_chars = tmp_left_chars
								right_chars = tmp_right_chars
								tsd_percent = round(100 * match_sum / totl_size, 1)
								
		if winning_tsd > -1:
			has_tsd = True
			
			
		return has_tsd, left_chars, right_chars, tsd_percent
		
		