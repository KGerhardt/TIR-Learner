import os
import warnings
import json

warnings.filterwarnings("ignore", category=UserWarning)  # mute keras warning

import torch                                                                    # noqa
import keras                                                                    # noqa

from .new_seq_reader import json_loader, bed_worker, json_structure

import multiprocessing

import numpy as np

PROGRAM_ROOT_DIR_PATH = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
CNN_MODEL_DIR_ABS_PATH = os.path.join(PROGRAM_ROOT_DIR_PATH, 'cnn0912', 'cnn0912.keras')

global feature_size
feature_size = 200
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
	cnn_seqs = []
	
	#These are tuples of strings as (tsd1, tir1, any_middle_sequence, tir2, tsd2)
	#All of these are in the FORWARD orientation, so tir2 needs RC for checking later
	
	for seq_tuple in bl.my_loaded_sequences:
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
		
		#Encode as integers
		arr = np.array([encoder[c] for c in cnn_sequence]).reshape(-1, 1)
		
		#Prepare as one-hot encoding for CNN
		hotboi = keras.utils.to_categorical(arr, num_classes=5)
		
		#Store the sequences for search
		cnn_seqs.append(hotboi)
	

	clean_json, final_gff3, final_fasta, keep_indices = None, None, None, None
	
	if len(cnn_seqs) > 0:
		#Final prep + model search
		cnn_seqs = np.stack(cnn_seqs)
				
		predicted_labels = model.predict(cnn_seqs, verbose = None)
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
			clean_json, final_gff3, final_fasta, keep_indices = bl.cnn_filter_json(passing_indices, 
																					numpy_classes, 
																					tsd_percentages,
																					tir_percentages,
																					module = 'Module4',
																					cnn_scores = retained_cnn_labels)
																					#cnn_scores = None)
			
		else:
			clean_json, final_gff3, final_fasta, keep_indices = None, None, None, None

	return bl.source, clean_json, final_gff3, final_fasta, keep_indices
	
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
		
		#We have four target final files, which can be produced directly in this step
		#TIR-Learner_FinalAnn.fa - passing CNN FASTA
		#TIR-Learner_FinalAnn.gff3 - passing CNN GFF3
		#TIR-Learner_FinalAnn_filter.fa - passing CNN FASTA overlaps resolved
		#TIR-Learner_FinalAnn_filter.gff3 - passing CNN GFF3 overlaps resolved
		
		#I am including 4 more files:
		#GRF output JSON.gz
		#TIRVish output JSON.gz
		#CNN filtered TIRVISH JSON.gz
		#CNN filtered GRF JSON.gz
		
		final_fa = os.path.join(wd, 'current_results', 'TIR-Learner_FinalAnn.fa')
		final_g3 = os.path.join(wd, 'current_results', 'TIR-Learner_FinalAnn.gff3')
		final_fa_filt = os.path.join(wd, 'current_results', 'TIR-Learner_FinalAnn_filter.fa')
		final_g3_filt = os.path.join(wd, 'current_results', 'TIR-Learner_FinalAnn_filter.gff3')
		
		final_tirvish = None
		final_grf = None
		
		o1 = open(final_fa, 'w')
		o2 = open(final_g3, 'w')
		o3 = open(final_fa_filt, 'w')
		o4 = open(final_g3_filt, 'w')
		
		loader = json_loader(working_dir = wd)
		
		if self.tirvish is not None:
			post_cnn_file = os.path.join(wd, 'checkpoints', 'post_CNN_TIRVish_json.txt')
			if not os.path.exists(post_cnn_file):
			
				loader.load_json_for_cnn(self.tirvish)
				total_sequences = sum(len(w[1]) for w in loader.workloads)
				num_args = len(loader.workloads)
				ct = 0
				percent_mod = int((num_args / 100)+0.5) if num_args > 100 else 1
				
				print(f'Beginning TIRVish TIR candidate CNN search...')
				print('')
				
				final_json = {}
				with multiprocessing.Pool(self.threads, initializer = cnn_init, initargs = ()) as pool:
					for src, ret_json, gff3, fasta, keeps in pool.imap_unordered(one_cnn, loader.workloads):
						#print(gff3)
						ct += 1
						if ret_json is not None:
							for seqid in keeps:
								#Add keep/remove overlap info
								ret_json[seqid]['sequence_kept_after_overlaps'] = keeps[seqid].astype(int).tolist()
								
								#Print fasta, gff, and filtered fasta, gff outputs
								for i in range(0, keeps[seqid].shape[0]):
									print(fasta[seqid][i], file = o1)
									print(gff3[seqid][i], file = o2)
									if keeps[seqid][i]:
										print(fasta[seqid][i], file = o3)
										print(gff3[seqid][i], file = o4)

							final_json[src] = ret_json
						if ct % percent_mod == 0:
							print(f'CNN search of TIRVish TIR candidates is {round(100*ct/num_args, 2)}% complete ({ct} of {num_args})')
						
				print('')
				print('Writing TIRVish JSON output...')
				with open(post_cnn_file, 'w', encoding = 'ascii') as out:
					json.dump(final_json, out, indent = 4)
				
				print(f'TIRVish TIR candidate CNN search complete!')
				print('')
			else:
				print('TIRVish CNN search already complete.')
				
			final_tirvish = post_cnn_file

		if self.grf is not None:
			loader.load_json_for_cnn(self.grf)
			num_args = len(loader.workloads)
			ct = 0
			percent_mod = int((num_args / 100)+0.5) if num_args > 100 else 1
			
			post_cnn_file = os.path.join(wd, 'checkpoints', 'post_CNN_GRF_json.txt')
			
			if not os.path.exists(post_cnn_file):
				print(f'Beginning GRF TIR candidate CNN search...')
				print('')
				
				final_json = {}
				with multiprocessing.Pool(self.threads, initializer = cnn_init, initargs = ()) as pool:
					for src, ret_json, gff3, fasta, keeps in pool.imap_unordered(one_cnn, loader.workloads):
						#print(gff3)
						if ret_json is not None:
							for seqid in keeps:
								#Add keep/remove overlap info
								ret_json[seqid]['sequence_kept_after_overlaps'] = keeps[seqid].astype(int).tolist()
								
								#Print fasta, gff, and filtered fasta, gff outputs
								for i in range(0, keeps[seqid].shape[0]):
									print(fasta[seqid][i], file = o1)
									print(gff3[seqid][i], file = o2)
									if keeps[seqid][i]:
										print(fasta[seqid][i], file = o3)
										print(gff3[seqid][i], file = o4)
							
							#Store up the final json for one output at the end
							final_json[src] = ret_json
						ct += 1
						if ct % percent_mod == 0:
							print(f'CNN search of GRF TIR candidates is {round(100*ct/num_args, 2)}% complete ({ct} of {num_args})')
				
				print('')
				print('Writing GRF JSON output...')
				with open(post_cnn_file, 'w', encoding = 'ascii') as out:
					json.dump(final_json, out, indent = 4)
					
				print(f'GRF TIR candidate CNN search complete!')
				print('')
			else:
				print('GRF CNN search already complete.')
			
			final_grf = post_cnn_file
		
		o1.close()
		o2.close()
		o3.close()
		o4.close()
		
		return final_fa, final_g3, final_fa_filt, final_g3_filt, final_tirvish, final_grf
		
		
		
		
			
			

'''
def run_by_chunks(chunk):
	df = pd.read_csv(chunk, sep = "\t", dtype = {'id':str, 'seq':str}, names = ("id", "seq"))

	#print("  Step 1/7: Getting sequence fragment for prediction")
	#df["seq_frag"] = df.swifter.progress_bar(False).apply(_get_sequence_fragment, axis=1)
	#just process in-flow rather than loading the data again later
	#df = df.drop(columns="seq")

	#df = _feature_encoding(df, False)
	
	df = _ml_predict(df, CNN_MODEL_DIR_ABS_PATH)
	
	if df is not None:
		df = _post_processing(df, False)
	
	os.remove(chunk)
	
	return df


def _get_sequence_fragment(x: pd.Series, feature_size: int = 200) -> str:
	#200 bp extended sequences
	seq = x["seq"]
	len_seq = len(seq)
	
	#Since these are the extended sequences, this would just be the flanks with no TIR at all
	if len_seq >= feature_size * 2:
		return seq[0:feature_size] + seq[-feature_size:]
	
	half_seqlen = int(len_seq / 2)
	
	#Front half of the sequence
	s1 = seq[0:half_seqlen]

	#Back half of the sequence
	s2 = seq[half_seqlen:]
	
	#Padding that extend s1 and s2 to feature size
	n1 = "N" * (feature_size - len(s1))
	n2 = "N" * (feature_size - len(s2))
	
	#I think these should be:
	#s1 = n1 + s1
	#s2 = 
	
	s1 += n1
	s2 = n2 + s2
	
	#Is this really correct? this will be s1NNNNNNNNs2NNNNNNNNNNN
	#Wouldn't we expect NNNNNNNS1S2NNNNNNNNN?
	return s1 + s2


def _feature_encoding(df_in: pd.DataFrame, flag_verbose: bool = True) -> pd.DataFrame:
	feature_int_encoder = LabelEncoder()
	voc = ["A", "C", "G", "T", "N"]
	num_classes = len(voc)
	feature_int_encoder.fit(voc)

	df = df_in.loc[:, ["id", "seq_frag"]].copy()
	#print("  Step 2/7: Label Encoding - Transforming non-numerical labels to numerical labels")
	df["int_enc"] = df.swifter.progress_bar(flag_verbose).apply(
		lambda x: np.array(feature_int_encoder.transform(list(x["seq_frag"]))).reshape(-1, 1), axis=1)
	df = df.drop(columns="seq_frag")

	#print("  Step 3/7: One-Hot Encoding - Converting class vectors to binary class matrices")
	df["feature"] = df.swifter.progress_bar(flag_verbose).apply(
		lambda x: keras.utils.to_categorical(x["int_enc"], num_classes=num_classes), axis=1)
	df = df.drop(columns="int_enc")

	return df



#This will get removed later

from typing import Set, Tuple, List, Dict, Iterable, Sequence, Union, Optional  # noqa

#def _ml_predict(df_in: pd.DataFrame, genome_file: str, path_to_model: str) -> Optional[pd.DataFrame]:
def _ml_predict(df_in: pd.DataFrame, path_to_model: str) -> Optional[pd.DataFrame]:
	model = keras.models.load_model(path_to_model)
	pre_feature = df_in["feature"].to_numpy()
	df = df_in.drop(columns="feature")

	if pre_feature.shape[0] == 0:
		#print("Info: " + genome_file + " has no candidate to be classified")
		return None

	l_class = ["DTA", "DTC", "DTH", "DTM", "DTT", "NonTIR"]
	target_int_encoder = LabelEncoder()
	target_int_encoder.fit(l_class)
	target_int_encoded = target_int_encoder.transform(l_class)
	d = dict(zip(target_int_encoded, l_class))

	#print("  Step 4/7: CNN prediction")
	predicted_labels = model.predict(np.stack(pre_feature), verbose = None)

	df["percent"] = pd.Series(predicted_labels.max(axis=-1))
	y_classes = predicted_labels.argmax(axis=-1)
	df["TIR_type"] = pd.Series([d[i] for i in y_classes])
	
	return df




def _post_processing(df_in: pd.DataFrame, flag_verbose: bool = True) -> pd.DataFrame:
	df = df_in.loc[:, ["id", "TIR_type"]]
	df = df[df["TIR_type"] != "NonTIR"].reset_index(drop=True)
	#print("  Step 5/7: Retrieving sequence ID")
	df["seqid"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: x["id"].split(":")[0], axis=1)
	#print("  Step 6/7: Retrieving sequence starting coordinate")
	df["sstart"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: int(x["id"].split(":")[1]), axis=1)
	#print("  Step 7/7: Retrieving sequence ending coordinate")
	df["send"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: int(x["id"].split(":")[2]), axis=1)
	df = df.loc[:, ["TIR_type", "id", "seqid", "sstart", "send"]]
	df = df.sort_values(["TIR_type", "seqid", "sstart", "send"], ignore_index=True)
		
	return df

def execute(TIRLearner_instance) -> pd.DataFrame:
	num_chunks = min([TIRLearner_instance.processors, len(TIRLearner_instance["base"])])

	ct = 0
	print(f"{len(TIRLearner_instance["base"])} partitions to process with CNN")
	pool = mp.Pool(num_chunks)
	with open(TIRLearner_instance.processed_de_novo_result_file_name_cnn, "w", newline='') as out:
		for result in pool.imap_unordered(run_by_chunks, TIRLearner_instance["base"]):
			ct += 1
			if result is not None:
				result = _get_start_end(TIRLearner_instance.genome_lengths, result)
				result.to_csv(out, sep = '\t', header = False, index = False)
			print(f'{ct} of {len(TIRLearner_instance["base"])} complete!')

	return TIRLearner_instance.processed_de_novo_result_file_name_cnn
	
'''


'''
dir_base = 'bigboi/'
tirvish = dir_base + 'checkpoints/TIRVish_json.txt'
grf = dir_base + 'checkpoints/GRF_json.txt'
work = dir_base + 'current_results/'

print('Running CNN')
mn = CNN_manager(tirvish, grf, work, threads)
mn.run()
'''