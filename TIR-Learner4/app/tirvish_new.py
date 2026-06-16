import os
import multiprocessing
import subprocess
import shutil

from .get_tans import tan_worker
from .new_tir_tsd import tsd_tir_checker
from .new_seq_reader import json_structure, dereplicate_json
from .output_compressor import write_formatted_json

#Rust TIRvish detection binary (drop-in for gt suffixerator + gt tirvish). Override via TIRVISH_RS_BIN env var.
#The conda package installs it as both `tirvish` and `tirvish_rs`; default to `tirvish`.
TIRVISH_RS_BIN = os.environ.get('TIRVISH_RS_BIN', 'tirvish')

def tirvish_init(output_dir, max_tir_len):
	global outdir
	outdir = output_dir
	global TIR_length
	TIR_length = max_tir_len

def one_tirvish(gen):
	chunk_base = os.path.basename(gen)

	#tirvish_rs --batch writes one flat TSV per chunk: <outdir>/<chunk_stem>.tirvish.tsv (no per-chunk dir)
	actual_output_file = os.path.join(outdir, f'{os.path.splitext(chunk_base)[0]}.tirvish.tsv')

	#Detection is done up-front by one tirvish_rs --batch call in TIRvish_manager; the TSV already exists.

	max_N_pct = 0.2
	max_ta_pct = 0.7

	#keep_sequences: retain each chunk's sequence in memory so each TIR arm can be re-sliced by coordinate.
	tan_check = tan_worker(gen, keep_sequences = True)
	aligner = tsd_tir_checker()
	json_record = json_structure(tan_check.seqlens, include_label = False)

	#tirvish_rs emits one row per hit carrying exactly the six (start, stop) coordinate pairs gt tirvish
	#writes per element -- the same data the old GFF parse collected into `next_result[0..5]`, but as plain
	#columns instead of six feature lines + regex. Columns (1-based, inclusive, local to the contig):
	#  0 seqid
	#  1,2   full element (repeat_region)            -> next_result[0]
	#  3,4   TSD1                                     -> next_result[1]
	#  5,6   body (terminal_inverted_repeat_element)  -> next_result[2]
	#  7,8   TIR1   (gt's (start,end)-first arm)      -> next_result[3]
	#  9,10  TIR2                                     -> next_result[4]
	#  11,12 TSD2                                     -> next_result[5]
	#  13    sim (tir_similarity; diagnostic, ignored here -- wfa_align re-derives similarity)
	#We rebuild next_result from the columns and then run the IDENTICAL calculations the GFF path ran, so
	#behaviour is byte-for-byte unchanged -- only the ingestion is simpler. seqid is the genomeSplitter
	#"{chrom};;{offset}" key, identical to the tan_check / seq_dict key (carried through verbatim).
	with open(actual_output_file) as tf:
		tf.readline()  # discard the column header
		for line in tf:
			line = line.rstrip('\n')
			if not line:
				continue
			f = line.split('\t')
			seqid = f[0]
			next_result = [(int(f[1]),  int(f[2])),    # [0] full element
						   (int(f[3]),  int(f[4])),    # [1] TSD1
						   (int(f[5]),  int(f[6])),    # [2] body (no TSD)
						   (int(f[7]),  int(f[8])),    # [3] TIR1
						   (int(f[9]),  int(f[10])),   # [4] TIR2
						   (int(f[11]), int(f[12]))]   # [5] TSD2

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

			#Check TA% and N% for the element body and each TIR arm, WITHOUT the TSD (as before)
			ok_seq = tan_check.check_acceptable_tans(seqid, no_tsd_start, no_tsd_stop, min_seqlen = 0, max_ta_pct = max_ta_pct, max_n_pct = max_N_pct)
			ok_tir1 = tan_check.check_acceptable_tans(seqid, tir_1_start, tir_1_stop, min_seqlen = 0, max_ta_pct = max_ta_pct, max_n_pct = max_N_pct)
			ok_tir2 = tan_check.check_acceptable_tans(seqid, tir_2_start, tir_2_stop, min_seqlen = 0, max_ta_pct = max_ta_pct, max_n_pct = max_N_pct)

			if ok_seq and ok_tir1 and ok_tir2:
				#TIRvish doesn't include the TSD in the element body, but full start/stop already carry it
				#(repeat_region), which is what we want downstream; so no TSD adjustment is needed here.

				#Check TIR length + similarity so we can pre-filter before the CNN
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

	os.remove(actual_output_file)

	json_record.sort_records()

	return json_record, gen

def TIRvish_manager(input_genome_files, original_genome_seqlen_dict, output_directory, checkpoint_directory, overlap_size, chunk_size, threads = 1, max_TIR_length = 5000):
	checkf = os.path.join(checkpoint_directory, 'TIRVish_json.txt')
	tirvish_json = os.path.join(output_directory, 'TIRVish_json.txt')

	#If the output exists, skip
	if not os.path.exists(checkf):
		args = input_genome_files

		num_args = len(args)
		print(f'There are {num_args} inputs to process with TIRVish')

		ct = 0
		#Give the user no more than 100 updates
		percent_mod = int((num_args / 100)+0.5) if num_args > 100 else 1

		combined_json = {}

		#--- Rust TIRvish detection: ONE work-stealing batch over all chunks (replaces per-file gt) ---
		#tirvish_rs --batch <outdir> reads the fragment paths from stdin and writes
		#<outdir>/<chunk_stem>.tirvish.tsv per fragment, parallel at the fragment level (1 frag/worker,
		#tail seed-steal). It runs directly on the FASTA, so no `gt suffixerator` index step is needed.
		listfile = os.path.join(output_directory, 'tirvish_batch_list.txt')
		with open(listfile, 'w') as lf:
			for g in args:
				lf.write(g + '\n')
		batch_cmd = [TIRVISH_RS_BIN, '--batch', output_directory, '--threads', str(threads),
			'-seed', '20', '-mintirlen', '10', '-maxtirlen', '1000',
			'-mintirdist', '10', '-maxtirdist', str(max_TIR_length),
			'-similar', '80', '-mintsd', '2', '-maxtsd', '11', '-vic', '13']
		with open(listfile) as lf:
			subprocess.run(batch_cmd, stdin=lf, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

		with multiprocessing.Pool(threads, initializer=tirvish_init, initargs=(output_directory, max_TIR_length)) as pool:
			for json_dict, genome_file in pool.imap_unordered(one_tirvish, args):
				if json_dict.has_records:
					combined_json[genome_file] = json_dict.json_record
				ct += 1
				if ct % percent_mod == 0:
					print(f'TIRVish search is {round(100*ct/num_args, 2)}% complete ({ct} of {num_args})')

				#Clean up: one_tirvish already removed the flat candidate TSV; no per-chunk dir to reap.

		combined_json = dereplicate_json(combined_json, overlap_size)

		write_formatted_json(combined_json, tirvish_json)

		#Checkpoint code here
		shutil.copy(tirvish_json, checkf)

	return checkf
