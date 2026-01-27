import os
import multiprocessing
import subprocess
import re
import shutil
import json

from .get_tans import tan_worker
from .new_tir_tsd import tsd_tir_checker
from .new_seq_reader import json_structure, dereplicate_json

#Regular expressions I use later on
genome_split_regex = re.compile(r'(.+);;(\d+)')
long_chunk_offset_regex = re.compile(r'long_chunk_(.+)_offset_(\d+).fasta')

def one_tirvish(args):
	genome = args[0]
	idx = args[1]
	
	#Build the TIRvish suffix array
	comm = f'gt suffixerator -db {genome} -indexname {idx} -tis -suf -lcp -des -ssp -sds -dna -mirrored'
	comm = comm.split()
	proc = subprocess.run(comm)
	
	#Execute the TIRvish search
	comm = ' '.join([f'gt tirvish -index {idx} -seed 20 -mintirlen 10 -maxtirlen 1000',
		   f'-mintirdist 10 -maxtirdist {5000} -similar 80 -mintsd 2 -maxtsd 11',
		   f'-vic 13 -seqids \"yes\"'])
		
	comm = comm.split()
	proc = subprocess.run(comm, capture_output = True, text = True)

	#Clean up the suffix array files
	for extension in ['des', 'esq', 'lcp', 'llv', 'md5', 'prj', 'sds', 'suf', 'ssp']:
		this_file = f'{idx}.{extension}'
		if os.path.exists(this_file):
			os.remove(this_file)
	
	max_N_pct = 0.2
	max_ta_pct = 0.7
	
	#Check TA and N percent OK using this
	tan_check = tan_worker(genome, keep_sequences = True)
	aligner = tsd_tir_checker()
	json_record = json_structure(tan_check.seqlens, include_label = False)
	
	next_result = []
	
	if len(proc.stdout) > 0:
		for line in proc.stdout.splitlines():
			if line.startswith('#'):
				if len(next_result) > 0:
					mat = re.match(genome_split_regex, seqid)
					mat = mat.groups()
					offset = int(mat[1])
					short_id = mat[0]
					
					tsd_1_start, tsd_1_stop = next_result[1]
					tsd_2_start, tsd_2_stop = next_result[5]
					
					tir_1_start, tir_1_stop = next_result[3]
					tir_2_start, tir_2_stop = next_result[4]
					
					tsd1_size = tsd_1_stop - tsd_1_start + 1
					tsd2_size = tsd_2_stop - tsd_2_start + 1
					tir1_size = tir_1_stop - tir_1_start + 1
					tir2_size = tir_2_stop - tir_2_start + 1
					
					full_element_start, full_element_stop = min(next_result[0]), max(next_result[0])
					no_tsd_start, no_tsd_stop = min(next_result[2]), max(next_result[2])
					
					#Skip conditions related to chunking boundaries, extension
					
					#If the sequence passes we adjust start/stop by TSD size, but we don't want to do that here 
					#because we need to check the sequence TA and N% without including the TSD, yet
					
					#The recovered sequence isn't purged by being too close to an edge
					
					#Check TA% and N% for the whole seqience
					ok_seq = tan_check.check_acceptable_tans(seqid, no_tsd_start, no_tsd_stop, min_seqlen = 0, max_ta_pct = max_ta_pct, max_n_pct = max_N_pct)
					ok_tir1 = tan_check.check_acceptable_tans(seqid, tir_1_start, tir_1_stop, min_seqlen = 0, max_ta_pct = max_ta_pct, max_n_pct = max_N_pct)
					ok_tir2 = tan_check.check_acceptable_tans(seqid, tir_2_start, tir_2_stop, min_seqlen = 0, max_ta_pct = max_ta_pct, max_n_pct = max_N_pct)
					
					if ok_seq and ok_tir1 and ok_tir2:
						
						#TIRvish doesn't include the TSD in the default start/end, but we want it for later purposes,
						#So we adjust start/stop coords to include it
						#seq_start = seq_start - tsd_size
						#seq_stop = seq_stop + tsd_size
						
						#Check TIR is acceptable length and similarity so that we don't have to post CNN and can filter out sequences
						left_tir_seq = tan_check.seq_dict[seqid][tir_1_start-1 : tir_1_stop]
						right_tir_seq = aligner.revcomp(tan_check.seq_dict[seqid][tir_2_start-1 : tir_2_stop])
						
						
						has_tir, l_rep_sz, r_rep_sz, r_start, q_start, pct = aligner.wfa_align(left_tir_seq, right_tir_seq, 
																						min_size = 10, min_similarity = 0.8)
						
						if has_tir:
							#Add JSON data
							json_record.add_record(seqid = seqid, 
													start = full_element_start, 
													stop = full_element_stop, 
													tsd1 = tsd1_size, 
													tsd2 = tsd2_size, 
													tir1 = tir1_size, 
													tir2 = tir2_size)
							

				next_result = []
			else:
				if 'tir_similarity=' in line:
					my_desc = line.strip()
				
				segs = line.strip().split('\t')
				seqid = segs[0]
				
				start, end = int(segs[3]), int(segs[4])
				#start, end = start + offset, end + offset
				
				typ = segs[2]
				
				next_result.append((start, end,))	
		
	json_record.sort_records()
	
	return json_record, genome

def TIRvish_manager(input_genome_files, original_genome_seqlen_dict, output_directory, checkpoint_directory, overlap_size, chunk_size, threads = 1):
	#global olap_size
	#olap_size = overlap_size
	global chunk_sz
	chunk_sz = chunk_size
	global original_seqlens
	original_seqlens = original_genome_seqlen_dict
	global outdir
	outdir = output_directory

	outf = os.path.join(output_directory, 'TIRVish_results.txt')
	checkf = os.path.join(checkpoint_directory, 'TIRVish_json.txt')
	tirvish_json = os.path.join(output_directory, 'TIRVish_json.txt')
	
	#If the output exists, skip
	if not os.path.exists(checkf):
		base_gens = [os.path.basename(g) for g in input_genome_files]
		tirvish_indices = [os.path.join(output_directory, f'{g}.tirvish_idx') for g in base_gens]
		
		args = list(zip(input_genome_files, tirvish_indices))

		num_args = len(args)
		print(f'There are {num_args} inputs to process with TIRVish')

		ct = 0
		#Prevent more than 100 updates
		percent_mod = int((num_args / 100)+0.5) if num_args > 100 else 1
		
		combined_json = {}
		
		#with open(outf, 'wb') as out:
		with multiprocessing.Pool(threads) as pool:
			#for seq_file, json_dict, genome_file in pool.imap_unordered(one_tirvish, args):
			for json_dict, genome_file in pool.imap_unordered(one_tirvish, args):
				if json_dict.has_records:
					combined_json[genome_file] = json_dict.json_record
				ct += 1
				if ct % percent_mod == 0:
					print(f'TIRVish search is {round(100*ct/num_args, 2)}% complete ({ct} of {num_args})')
		
		#with open(tirvish_json+'_prefilter.txt', 'w', encoding = 'ascii') as out:
		#	json.dump(combined_json, out, indent = 4)
		
		combined_json = dereplicate_json(combined_json, overlap_size)
		
		with open(tirvish_json, 'w', encoding = 'ascii') as out:
			json.dump(combined_json, out, indent = 4)
		
		#Checkpoint code here
		shutil.copy(tirvish_json, checkf)
		
	return checkf