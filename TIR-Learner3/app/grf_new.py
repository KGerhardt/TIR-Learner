import sys
import os
import multiprocessing
import pyfastx
import subprocess

import re

import json
import shutil

from .get_tans import tan_worker
from .new_tir_tsd import tsd_tir_checker
from .new_seq_reader import json_structure, dereplicate_json

#Regular expressions used here and there
genome_split_regex = re.compile(r'(.+);;(\d+)')
grp_regex = re.compile(r'_offset_(\d+).fasta')
#>NC_069499.1;;693555512:491296:491325:8m1M1m:CG
grfmite_regex = re.compile(r'(.+);;(\d+):(\d+):(\d+):(.+):(.+)')
#long_chunk_NC_069499.1_offset_0.fasta
long_chunk_offset_regex = re.compile(r'long_chunk_(.+)_offset_(\d+).fasta')

#This one is pretty simple
cig_parse_regex = re.compile(r'(\d+)')

def parse_cig(GRF_cigar_string):
	#The cigar string-like component represents the TIR sequence. From the github:
	#The call for GRF specifically has --max_index=0, which prohibits any TIR indel
	#therefore only m/M are needed to be wrangled, each of which represents a sequence match and a TIR extension
	#politely, this means TIRs are symmetrical
	'''
	For TIR_pairing:
		"m" means matches in base pairing (A and T; C and G);
		"M" means mismatches in base pairing;
		"I" means insertions in base pairing;
		"D" means deletions in base pairing;
	'''
	
	tir_size = sum([int(i) for i in re.findall(cig_parse_regex, GRF_cigar_string)])
	
	return tir_size
	
def one_GRF(gen):
	#gen = args[0]
	
	chunk_base = os.path.basename(gen)
	
	out = os.path.join(outdir, f'{chunk_base}_grfmite')
	
	#out = args[1]
	actual_output_file = os.path.join(out, 'candidate.fasta')
	#partial_json = os.path.join(out, 'grf_partial_json.txt')
	
	comm = ' '.join([f'grf-main -i {gen} -o {out} -c 1 -t 1 -p 20 --min_space 10 --max_space {TIR_length}',
	f'--max_indel 0 --min_tr 10 --min_spacer_len 10 --max_spacer_len {TIR_length}'
	])
	
	#print(comm)
	comm = comm.split()
	
	#GRF command is very slow
	subprocess.run(comm, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
	
	min_seqlen = 50
	max_N_pct = 0.2
	max_ta_pct = 0.7
	
	#Local is_TA and is_N values
	#is_TA, is_N = get_tans(gen)
	tan_check = tan_worker(gen)
	aligner = tsd_tir_checker()
	json_record = json_structure(tan_check.seqlens, include_label = False)
	
	if 'short_chunk_' in gen:
		offset = 0
	else:
		#Need a check for "is the last sequence" which stops long checks from truncating too much
		mat = re.match(long_chunk_offset_regex, os.path.basename(gen)).groups()
		#It's always one sequence/file for a long chunk
		parent_seqid = mat[0]
		offset = int(mat[1])
	
	'''
	#Truncation checking goes here
	sequence_lengths = tan_check.seqlens
	
	needs_start_truncation = False
	needs_end_truncation = False
	
	start_truncation_cutoff = None
	end_truncation_cutoff = None
		
	#Short chunks always have 0 offset, do not require any truncation so the default needs start/end trunc = false prevents it later
	if 'short_chunk_' in gen:
		offset = 0
	else:
		#Need a check for "is the last sequence" which stops long checks from truncating too much
		mat = re.match(long_chunk_offset_regex, os.path.basename(gen)).groups()
		#It's always one sequence/file for a long chunk
		parent_seqid = mat[0]
		offset = int(mat[1])
		
		this_seqlen = sequence_lengths[f'{parent_seqid};;{offset}']
		original_genome_size = original_seqlens[parent_seqid]
		
		#We can't start truncate a sequence beginning at zero; a TIR found in the first 200bp of a chromosome is just wherever and whatever it is
		if offset > 0:
			needs_start_truncation = True
			start_truncation_cutoff = ext_size
			
		#Likewise, we can't truncate a TIR at the end of a sequence
		if offset + chunk_sz < original_genome_size:
			needs_end_truncation = True
			end_truncation_cutoff = this_seqlen - ext_size
	
		#For non-edge cases, we deliberately remove matches found in the first/last [EXTENSION_SIZE] bp of the chunk; 
		#these are not extensible within that chunk, but will always be found AND be extensible in another chunk thanks to the overlap
	'''

	#with open(cleaned_output_file, 'w') as outf:
	for seqid, sequence in pyfastx.Fasta(actual_output_file, build_index = False):
		first_4 = sequence[0:4]
		
		slo = len(sequence)
		
		match = re.match(grfmite_regex, seqid)
		grps = match.groups()
		
		#GRF TSDs are always symmetrical, which makes life easy on the JSON stuff
		tsd = grps[5]
		tsd_size = len(tsd)
		
		short_id = grps[0]
		
		#cig is the TIR element that needs parsed into tir1 start + stop and tir2 start + stop
		cig = grps[4]

		#Offset we can get from the file as above, not needed here it's a guaranteed constnat
		#offset = int(grps[1])
		
		#TA key will also be the seqid key inside the JSON
		ta_key = f'{short_id};;{offset}'
		
		#Note: GRF coordinates are not inclusive of the TSD; start-stop includes TIRs but not TSD
		
		#GRF is 1-indexed, but so is the tan_checker
		#The local_start -= tsd_size : local_stop += tsd_size is the full element with TSDs
		local_start, local_stop = int(grps[2]), int(grps[3])
		
		#these are the with-TSDs boundaries
		full_element_start, full_element_stop = local_start - tsd_size, local_stop + tsd_size
		
		#Uses prefix sum arrays to check if TA and N content, sequence length are OK
		#TA + N indices are intentionally 1-indexed to match the formatting of GRF
		
		#We use local start/stop here because we don't want to check the TSD TA/N percentage
		passing = tan_check.check_acceptable_tans(ta_key, local_start, local_stop, min_seqlen = min_seqlen, max_ta_pct = max_ta_pct, max_n_pct = max_N_pct)
		
		if passing:				
			#Passing sequences still need to check TSD and maybe first 4
			if tsd_size > 6 or tsd == "TAA" or tsd == "TTA" or tsd == "TA" or first_4 == "CACT" or first_4 == "GTGA":
				
				#TIRs are also symmetrical in GRF
				tir_size = parse_cig(cig)
				
				#Check TIR is acceptable length and similarity so that we don't have to post CNN
				left_tir_seq = sequence[0:tir_size]
				right_tir_seq = aligner.revcomp(sequence[-tir_size:])
				
				has_tir, l_rep_sz, r_rep_sz, r_start, q_start, pct = aligner.wfa_align(left_tir_seq, right_tir_seq, 
																				min_size = 10, min_similarity = 0.8)
				if has_tir:
					json_record.add_record(seqid = ta_key, 
											start = full_element_start, 
											stop = full_element_stop, 
											tsd1 = tsd_size, 
											tsd2 = tsd_size, 
											tir1 = tir_size, 
											tir2 = tir_size)
	
	
	
	os.remove(actual_output_file)

	json_record.sort_records()
	
	return json_record, gen, out
	
def GRF_manager(input_genome_files, original_genome_seqlen_dict, output_directory, checkpoint_directory, overlap_size, chunk_size, threads = 1, max_TIR_length = 5000):
	#global olap_size
	#olap_size = overlap_size
	global chunk_sz
	chunk_sz = chunk_size
	global original_seqlens
	original_seqlens = original_genome_seqlen_dict
	global outdir
	outdir = output_directory
	global TIR_length
	TIR_length = max_TIR_length

	#outf = os.path.join(output_directory, 'GRF_results.txt')
	checkf = os.path.join(checkpoint_directory, 'GRF_json.txt')
	grf_json = os.path.join(output_directory, 'GRF_json.txt')
	
	#If the output exists, skip
	if not os.path.exists(checkf):
		base_gens = [os.path.basename(g) for g in input_genome_files]
		#GRF_outputs = [os.path.join(output_directory, f'{g}.grfmite') for g in base_gens]
		#GRF_outputs = [os.path.join(output_directory, f'{g}.grfmite.candidate.fa') for g in base_gens]
		
		args = input_genome_files
		
		#args = list(zip(input_genome_files, GRF_outputs))
		#Here for testing
		#front = args[0:10]
		#back = args[-10:]
		#args = [*front, *back]
		#args = [args[0]]
		
		num_args = len(args)
		print(f'There are {num_args} inputs to process with GRF')

		ct = 0
		#Give the user no more than 100 updates
		percent_mod = int((num_args / 100)+0.5) if num_args > 100 else 1
		
		combined_json = {}
		#with open(outf, 'wb') as outfile:
		with multiprocessing.Pool(threads) as pool:
		#with multiprocessing.Pool(1) as pool:
			for json_dict, genome_file, output_directory in pool.imap_unordered(one_GRF, args):				
				#Convert out of numpy for json write
				if json_dict.has_records:
					combined_json[genome_file] = json_dict.json_record
				
				ct += 1
				if ct % percent_mod == 0:
					print(f'GRF search is {round(100*ct/num_args, 2)}% complete ({ct} of {num_args})')
					
				#Clean up
				#os.remove(output_tsv)
				os.rmdir(output_directory)
		
		combined_json = dereplicate_json(combined_json, overlap_size)
		
		with open(grf_json, 'w', encoding = 'ascii') as out:
			json.dump(combined_json, out, indent = 4)
			
		#Checkpoint code here
		shutil.copy(grf_json, checkf)	
	
	return checkf