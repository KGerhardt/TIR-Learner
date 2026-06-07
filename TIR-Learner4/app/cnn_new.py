import os
import warnings
import json

warnings.filterwarnings("ignore", category=UserWarning)  # mute keras warning

import torch                                                                    # noqa
import keras                                                                    # noqa

from .new_seq_reader import json_loader, bed_worker, json_structure
from .output_compressor import write_formatted_json

import multiprocessing

import numpy as np

PROGRAM_ROOT_DIR_PATH = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
CNN_MODEL_DIR_ABS_PATH = os.path.join(PROGRAM_ROOT_DIR_PATH, 'cnn0912', 'cnn0912.keras')

global feature_size
feature_size = 200

#Max candidates handed to one CNN worker batch. Caps per-worker RAM and gives
#real load balance (a 400k-candidate chunk becomes ~20 even batches, not one straggler).
CNN_BATCH_SIZE = 20000
global path_to_model

path_to_model = CNN_MODEL_DIR_ABS_PATH

#remove pandas later
#import pandas as pd
#The ML modules take fuckin forever to load

from .new_tir_tsd import tsd_tir_checker

#Load repeatedly used resources just once and reuse for each chunk
def cnn_init():
	global model
	model = keras.models.load_model(path_to_model)
	global repeat_checker
	repeat_checker = tsd_tir_checker()
	global encoder
	#Has to be alphabetical, but this is equivalent to the feature encoder from sklearn
	encoder = {'A':0, 'C':1, 'G':2, 'N':3, 'T':4}	
	global base_lut
	#256-entry byte -> nucleotide-class LUT built from `encoder`; unknown bytes fall
	#back to N (the symbol used for padding) so vectorized encoding can't KeyError.
	base_lut = np.full(256, encoder['N'], dtype = np.uint8)
	for _ch, _idx in encoder.items():
		base_lut[ord(_ch)] = _idx
	global onehot5
	#Identity rows give one-hot encoding by fancy-indexing; float32 matches model input.
	onehot5 = np.eye(len(encoder), dtype = np.float32)
	global l_class
	#Convert the l_class to a numpy array which enables easy numpy array slicing
	l_class = np.array(["DTA", "DTC", "DTH", "DTM", "DTT", "NonTIR"])

#The main function of this code
'''
(1) load candidate TIRs as tuples of strings encoding (tsd1, tir1, any_middle_sequence, tir2, tsd2)
(2) Collect the CNN search target from each - this is (up to) first and last 200 bp of tir1 + any_middle_sequence + tir2, OR
	pad first half of CNN seach target + (N * size) + second half of CNN search target to hit 400 bp total
(3) Format sequences for CNN search, execute CNN search	to get TIR class for each candidate or ID as NonTIR
(4) For all TIR candidates, check appropriate TIR start + end, TSD start + end patterns, TSD conservation% 
	(TA and N% checks, TIR conservation% check happened earlier in GRF / TIRVish filtering)
(5) Filter incoming JSONs to reflect only passing TIR sequences
(6) Format passing TIR sequences as FASTA, GFF, etc. as needed


'''
def one_cnn(workload):
	bl = bed_worker(workload)
	bl.load_refgen()
	bl.convert_json_to_sequences_for_cnn()
		
	#Select the first and last 200 bp and encode these for CNN model search
	n_candidates = len(bl.my_loaded_sequences)
	#Preallocate the (n, 400) integer-class matrix; each row is filled at C speed below.
	idx_mat = np.empty((n_candidates, 2 * feature_size), dtype = np.uint8)

	#These are tuples of strings as (tsd1, tir1, any_middle_sequence, tir2, tsd2)
	#All of these are in the FORWARD orientation, so tir2 needs RC for checking later

	for i, seq_tuple in enumerate(bl.my_loaded_sequences):
		#This is the sequence excluding TSDs
		my_tir_seqeuence = ''.join(seq_tuple[1:4])
		this_seqlen = len(my_tir_seqeuence)

		#Check if the seuqence needs padded
		if this_seqlen > 2 * feature_size:
			cnn_sequence = my_tir_seqeuence[0:feature_size] + my_tir_seqeuence[-feature_size:]
		else:
			#Ceiling function half sequence length
			left_size = int(this_seqlen / 2)
			right_size = this_seqlen - left_size
			
			N_size = 2 * feature_size - this_seqlen
			
			N_pad = 'N' * N_size
			
			#Pad with N's in the middle of what's there
			cnn_sequence = my_tir_seqeuence[0:left_size] + N_pad + my_tir_seqeuence[-right_size:]			
		
		#Map characters -> integer classes at C speed via the byte LUT. latin-1 is a
		#1:1 byte decode (identical to ascii for ACGNT); any non-ACGNT byte routes to
		#the LUT default (N) instead of raising, as the old encoder[c] lookup would.
		idx_mat[i] = base_lut[np.frombuffer(cnn_sequence.encode('latin-1'), dtype = np.uint8)]

	clean_json, final_gff3, final_fasta = None, None, None

	if n_candidates > 0:
		#One-hot the whole batch at once -> (n, 400, 5) float32. The planner caps a
		#batch at a modest candidate count, so this array stays small (a worker never
		#sees a 400k-candidate chunk anymore); batch_size below keeps predict's own
		#working set bounded and is ~2.7x faster than the keras default of 32.
		cnn_seqs = onehot5[idx_mat]

		predicted_labels = model.predict(cnn_seqs, verbose = None, batch_size = 256)
		#Free space
		cnn_seqs = None
		
		#We don't actually use the max per row data anywhere, so don't even bother
		#Pure numpy equivalents of the percent and class type selections
		#max_per_row = np.max(predicted_labels, axis = 1)
		
		#Select the CNN's assigned label for each sequence based on which probability is highest
		numpy_classes = np.argmax(predicted_labels, axis = 1)
		
		not_non_tirs = np.where(numpy_classes < 5)[0]
		
		passing_indices = []
		tir_percentages = []
		tsd_percentages = []
		
		for check_index, tir_type, in zip(not_non_tirs, l_class[numpy_classes[not_non_tirs]]):
			this_sequence = bl.my_loaded_sequences[check_index]
			
			tsd_1 = this_sequence[0]
			tir_1 = this_sequence[1]
			#This one gets rev-comp'd
			tir_2 = repeat_checker.revcomp(this_sequence[3])
			tsd_2 = this_sequence[4]
			
			#Check if either TIR sequence begins with the correct type of sequence
			ok_tir_conservation = repeat_checker.check_tir_conservation(tir_type, tir_1, tir_2)
			if ok_tir_conservation:
				#Check if the TIR has the correct type and acceptable size, similarity of TSD sequences
				has_tsd, left_tsd_size, right_tsd_size, tsd_percent = repeat_checker.check_tsd(tsd_1, 
																				tsd_2, 
																				tir_type = tir_type,
																				min_similarity = 0.8)
				#If so, collect;
				if has_tsd:
					#This re-alignment is strictly not necessary and a perfected version of this program 
					#would simply keep track of the percent through GRF and TIRVish by restructuring those JSON outputs;
					#I do not believe it's worthwhile currently, this is really only a few sequences and the process is fast 
					#and we already have the data loaded
					has_tir, l_rep_sz, r_rep_sz, r_start, q_start, pct = repeat_checker.wfa_align(tir_1, tir_2, 
																		min_size = 10, min_similarity = 0.8)
					
					tir_percentages.append(pct)
					tsd_percentages.append(tsd_percent)
					passing_indices.append(check_index)
					
		not_non_tirs = None
		passing_indices = np.array(passing_indices)
		
		if passing_indices.shape[0] > 0:
			numpy_classes = l_class[numpy_classes[passing_indices]].tolist()
			#retained_cnn_labels = predicted_labels[passing_indices, :].tolist()
			retained_cnn_labels = (np.round(predicted_labels[passing_indices, :], decimals = 4) * 10000).astype(np.int32).tolist()
			#for i in range(0, len(retained_cnn_labels)):
			#	retained_cnn_labels[i] = [ '%.5f' % lab for lab in retained_cnn_labels[i]]
				
			
			#Convert to a sequence-recovery JSON from the partial files and a 
			#dict of the final seqids and sequences with corrected positional indices
			#Overlap resolution is intentionally NOT done here: a batch holds only a
			#slice of a seqid's candidates, so overlaps are resolved in main once all
			#of a seqid's passers are merged (CNN_manager.run -> resolve_overlaps).
			clean_json, final_gff3, final_fasta = bl.format_cnn_passers(passing_indices,
																					numpy_classes,
																					tsd_percentages,
																					tir_percentages,
																					module = 'Module4',
																					cnn_scores = retained_cnn_labels)
																					#cnn_scores = None)

		else:
			clean_json, final_gff3, final_fasta = None, None, None

	return bl.source, clean_json, final_gff3, final_fasta
	
class CNN_manager:
	def __init__(self, tirvish, grf, working_dir = '.', threads = 1):
		global wd
		wd = working_dir
		
		self.tirvish = tirvish
		
		#Automatic no-homologs file checking for post-BLAST filtered module 1-3
		if self.tirvish is not None:
			homolog_file = os.path.join(wd, 'checkpoints', 'TIRVish_json_no_homologs.txt')
			#Check for no-homologs file automatically
			if os.path.exists(homolog_file):
				self.tirvish = homolog_file
		
		self.grf = grf
		#Automatic no-homologs file checking for post-BLAST filtered module 1-3
		if self.grf is not None:
			homolog_file = os.path.join(wd, 'checkpoints', 'GRF_json_no_homologs.txt')
			#Check for no-homologs file automatically
			if os.path.exists(homolog_file):
				self.grf = homolog_file
		
		self.threads = threads
		
		#self.run()
		
	def run(self):
		#Four target outputs, produced directly here:
		#  TIR-Learner_FinalAnn.fa / .gff3          - all passing-CNN candidates
		#  TIR-Learner_FinalAnn_filter.fa / .gff3   - passers after overlap resolution
		#plus the post-CNN JSON checkpoints (TIRVish, GRF).
		final_fa      = os.path.join(wd, 'current_results', 'TIR-Learner_FinalAnn.fa')
		final_g3      = os.path.join(wd, 'current_results', 'TIR-Learner_FinalAnn.gff3')
		final_fa_filt = os.path.join(wd, 'current_results', 'TIR-Learner_FinalAnn_filter.fa')
		final_g3_filt = os.path.join(wd, 'current_results', 'TIR-Learner_FinalAnn_filter.gff3')

		final_tirvish = None
		final_grf = None

		with open(final_fa, 'w') as o1, open(final_g3, 'w') as o2, open(final_fa_filt, 'w') as o3, open(final_g3_filt, 'w') as o4:
			if self.tirvish is not None:
				final_tirvish = os.path.join(wd, 'checkpoints', 'post_CNN_TIRVish_json.txt')
				self._cnn_search_branch(self.tirvish, final_tirvish, 'TIRVish', o1, o2, o3, o4)
			if self.grf is not None:
				final_grf = os.path.join(wd, 'checkpoints', 'post_CNN_GRF_json.txt')
				self._cnn_search_branch(self.grf, final_grf, 'GRF', o1, o2, o3, o4)

		return final_fa, final_g3, final_fa_filt, final_g3_filt, final_tirvish, final_grf

	def _cnn_search_branch(self, json_source, post_cnn_file, label, o_fa, o_g3, o_fa_filt, o_g3_filt):
		"""Run the CNN stage for one candidate source (TIRVish or GRF).

		Planner(main) -> CNN workers(batches) -> merge per seqid(main) -> overlap(main) -> emit.
		Workers run CNN + TSD/TIR checks on bounded candidate batches and return only passers
		with formatted GFF/FASTA. Overlap resolution is deferred to here so a seqid's passers
		- which the planner may split across batches - are resolved as a whole. Overlap consumes
		only coordinates and is order-dependent, so main re-sorts each seqid's merged passers with
		the exact upstream key (json_structure.sort_records) before resolving; this makes the
		result independent of batch arrival order, so workers can return out of order.
		"""
		if os.path.exists(post_cnn_file):
			print(f'{label} CNN search already complete.')
			return

		loader = json_loader(working_dir = wd)
		loader.load_json_for_cnn(json_source)
		batches = list(loader.build_cnn_batches(CNN_BATCH_SIZE))
		#Candidates are copied into batches; drop the full source JSON to cap main RAM.
		loader.json_data = None
		loader.workloads = None

		num_args = len(batches)
		ct = 0
		percent_mod = int((num_args / 100) + 0.5) if num_args > 100 else 1
		print(f'Beginning {label} TIR candidate CNN search...')
		print('')

		#Accumulate passers per source-chunk -> seqid (a seqid may span several batches).
		final_json = {}   # src -> seqid -> post-CNN record (per-candidate arrays + scalars)
		gff_acc    = {}   # src -> seqid -> [gff lines]     (aligned to the arrays)
		fasta_acc  = {}   # src -> seqid -> [fasta records] (aligned to the arrays)

		with multiprocessing.Pool(self.threads, initializer = cnn_init, initargs = ()) as pool:
			for src, clean_json, gff3, fasta in pool.imap_unordered(one_cnn, batches):
				ct += 1
				if clean_json is not None:
					fj = final_json.setdefault(src, {})
					ga = gff_acc.setdefault(src, {})
					fa = fasta_acc.setdefault(src, {})
					for seqid, block in clean_json.items():
						if seqid not in fj:
							fj[seqid] = block
							ga[seqid] = gff3[seqid]
							fa[seqid] = fasta[seqid]
						else:
							dst = fj[seqid]
							for k, v in block.items():
								#Extend per-candidate arrays; per-seqid scalars are already set.
								if isinstance(v, list):
									dst[k].extend(v)
							ga[seqid].extend(gff3[seqid])
							fa[seqid].extend(fasta[seqid])
				if ct % percent_mod == 0:
					print(f'CNN search of {label} TIR candidates is {round(100*ct/num_args, 2)}% complete ({ct} of {num_args})')

		print('')
		print(f'Resolving overlaps and writing {label} output...')
		for src in final_json:
			for seqid in final_json[src]:
				block = final_json[src][seqid]
				#Re-sort the merged passers into the exact upstream order before resolving
				#(matches json_structure.sort_records: primary start+tsd1 asc, secondary
				#stop+tsd2 desc). Overlap is order-dependent; this removes any dependence on
				#which order the workers happened to return batches in.
				starts    = np.array(block['seq_start_incl_tsd']) + np.array(block['tsd1_size'])
				neg_stops = -1 * np.array(block['seq_stop_incl_tsd']) - np.array(block['tsd2_size'])
				order = np.lexsort((neg_stops, starts)).tolist()
				for k, v in block.items():
					if isinstance(v, list):
						block[k] = [v[i] for i in order]
				gffs = [gff_acc[src][seqid][i] for i in order]
				fas  = [fasta_acc[src][seqid][i] for i in order]
				keep = bed_worker.resolve_overlaps({seqid: block})[seqid]
				block['sequence_kept_after_overlaps'] = keep.astype(int).tolist()
				for i in range(keep.shape[0]):
					print(fas[i], file = o_fa)
					print(gffs[i], file = o_g3)
					if keep[i]:
						print(fas[i], file = o_fa_filt)
						print(gffs[i], file = o_g3_filt)

		write_formatted_json(final_json, post_cnn_file)
		print(f'{label} TIR candidate CNN search complete!')
		print('')