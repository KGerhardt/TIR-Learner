#!/usr/app/env python3
# -*- coding: utf-8 -*-

from shared import *

from run_TIRvish import process_fasta, retrieve_split_sequence_offset

REGEX_PATTERN_PART_ORIGINAL_SEQ_ID = r"(.+)"
REGEX_PATTERN_PART_SEGMENT_POSITION_X = r"([\w.]+"
REGEX_PATTERN_PART_SEGMENT_POSITION_Y = r"\d+)"
SPLIT_SEQ_ID_REGEX_PATTERN = ('^' +
							  MP_SPLIT_SEQ_ID_FORMAT_STR.format(REGEX_PATTERN_PART_ORIGINAL_SEQ_ID,
																REGEX_PATTERN_PART_SEGMENT_POSITION_X,
																REGEX_PATTERN_PART_SEGMENT_POSITION_Y) +
							  r":(\d+):(\d+):(\w+):(\w+)$")
# SPLIT_SEQ_ID_PATTERN = r"^(\w+)_split_([\w.]+of\d+):(\d+):(\d+):(\w+):(\w+)$"

def _find_digits_sum(string: str) -> int:
	pattern = r"(\d+)"
	l = re.findall(pattern, string)
	return sum([int(i) for i in l])

def processors_allocation(processors: int) -> Tuple[int, int, int]:
	if processors <= 16:
		num_process = int(math.sqrt(processors))
	else:
		# Divide by sqrt(3) ≈ 1.732 for larger processor counts to favor more threads per process
		num_process = int(math.sqrt(processors/math.sqrt(3)))
	num_thread_per_process = processors // num_process
	total_num_threads = num_process * num_thread_per_process
	return num_process, num_thread_per_process, total_num_threads


def retrieve_unsplit_seq_id(match: re.Match, split_seq_len: int, overlap_seq_len: int) -> str:
	seqid, segment_position, start, end, TIR_pattern, TSD = match.groups()
	offset = retrieve_split_sequence_offset(segment_position, split_seq_len, overlap_seq_len)
	return f"{seqid}:{int(start) + offset}:{int(end) + offset}:{TIR_pattern}:{TSD}"

def _GRF(genome_file: str, genome_name: str, TIR_length: int, processors: int, GRF_path: str):
	GRF_bin_path = os.path.join(GRF_path, "grf-main")
	GRF_result_dir_name = f"{genome_name}_GRFmite"
	grf = (f"\"{GRF_bin_path}\" -i \"{genome_file}\" -o \"{GRF_result_dir_name}\" -c 1 -t {processors} -p 20 "
		   f"--min_space 10 --max_space {TIR_length} --max_indel 0 --min_tr 10 "
		   f"--min_spacer_len 10 --max_spacer_len {TIR_length}")
	shell_filter = r" | grep -vE 'start:|end:|print:|calculate|^$'"
	subprocess.Popen(grf + shell_filter, shell=True).wait()

def _GRF_mp(genome_file_path: str, genome_name: str, TIR_length: int, processors: int, GRF_path: str):
	GRF_working_dir = os.path.dirname(genome_file_path)
	genome_file_name = os.path.basename(genome_file_path)
	os.chdir(GRF_working_dir)
	_GRF(genome_file_name, genome_name, TIR_length, processors, GRF_path)
	os.chdir("../")

def _run_GRF_native(genome_file: str, genome_name: str, TIR_length: int,
					processors: int, flag_debug: bool, GRF_path: str) -> Optional[pd.DataFrame]:
	print("  Step 1/2: Executing GRF in native mode\n")
	_GRF(genome_file, genome_name, TIR_length, processors, GRF_path)
	print("  Step 2/2: Getting GRF result")
	return _get_GRF_result_df_native(genome_name, flag_debug)


def _run_GRF_py_para(genome_file: str, genome_name: str, TIR_length: int, processors: int,
					 flag_debug: bool, GRF_path: str, fasta_files_path_list: List[str], output_file: str) -> pd.DataFrame:
	os.makedirs(f"{SPLIT_FASTA_TAG}_mp", exist_ok=True)
	os.chdir(f"./{SPLIT_FASTA_TAG}_mp")

	print("  Step 1/3: Checking processed FASTA files")
	if (not fasta_files_path_list or
			any(not os.path.exists(f) or os.path.getsize(f) == 0 for f in fasta_files_path_list)):
		print("	Processed FASTA files not found / invalid. Re-process FASTA files.")
		fasta_files_path_list.extend(process_fasta(genome_file, MP_SPLIT_SEQ_LEN, MP_OVERLAP_SEQ_LEN))

	

	print("  Step 2/3: Executing GRF with python multiprocessing")
	num_process, num_thread_per_process, _ = processors_allocation(processors)

	mp_args_list = [(file_path, genome_name, TIR_length, num_thread_per_process, GRF_path) for
					file_path in fasta_files_path_list]
	
	print()
	with mp.Pool(num_process) as pool:
		pool.starmap(_GRF_mp, mp_args_list)
	print()

	print("  Step 3/3: Getting GRF result")
	
	df = _get_GRF_result_df_para(fasta_files_path_list, genome_name, processors, flag_debug,
									   MP_SPLIT_SEQ_LEN, MP_OVERLAP_SEQ_LEN, output_file)
	
	return df


'''
'''
#def _get_single_GRF_result_df_para(file_path: str, flag_debug: bool,
#								   split_seq_len: int, overlap_seq_len: int) -> Optional[pd.DataFrame]:
def _get_single_GRF_result_df_para(args) -> Optional[pd.DataFrame]:
	file_path, flag_debug, split_seq_len, overlap_seq_len = args
	seen_records = set()
	try:
		min_seqlen = 50
		max_N_pct = 0.2
		max_ta_pct = 0.7
		
		df_data_dict = {}
		for rec in SeqIO.parse(os.path.join(file_path, "candidate.fasta"), "fasta"):
			caps = str(rec.seq)
			caps = caps.upper()
			seqlen = len(caps)
			if seqlen > min_seqlen:
				cts = Counter(caps)
				if (cts['N'] / seqlen) < max_N_pct and ((cts['A']+cts['T']) / seqlen) < max_ta_pct:
					seqname = rec.id
					
					if MP_SPLIT_SEQ_TAG in seqname:
						match = re.match(SPLIT_SEQ_ID_REGEX_PATTERN, seqname)
						seqname = retrieve_unsplit_seq_id(match, split_seq_len, overlap_seq_len)
					
					tir_len = _find_digits_sum(seqname.split(":")[-2])
					tir = caps[0:tir_len]
					cts = Counter(tir)
					if (cts['N'] / tir_len) < max_N_pct and ((cts['A']+cts['T']) / tir_len) < max_ta_pct:
						tsd = seqname.split(":")[-1]
						first_4 = caps[0:4]
						if len(tsd) > 6 or tsd == "TAA" or tsd == "TTA" or tsd == "TA" or first_4 == "CACT" or first_4 == "GTGA":
							seen_records.add(seqname)
							df_data_dict[seqname] = caps
		
	except FileNotFoundError:
		df_data_dict = None
	if not flag_debug:
		subprocess.Popen(["rm", "-rf", file_path])
	return df_data_dict, seen_records
	
	
def _get_GRF_result_df_para(fasta_files_path_list: List[str], genome_name: str, processors: int, flag_debug: bool,
							split_seq_len: int, overlap_seq_len: int, output_file: str) -> Optional[pd.DataFrame]:
	GRF_result_dir_name = f"{genome_name}_GRFmite"
	GRF_result_dir_list = [os.path.join(os.path.dirname(file), GRF_result_dir_name) for file in
						   fasta_files_path_list]

	mp_args_list = [(file_path, flag_debug, split_seq_len, overlap_seq_len) for file_path in GRF_result_dir_list]

	df_list = []
	ct = 0
	argnum = len(mp_args_list)
	total_seen_records = set()
	if argnum > 100:
		print(f"{argnum} GRF results to process.")
	with open(output_file, "w") as out:
		with mp.Pool(processors) as pool:
			for df, seen_records in pool.imap_unordered(_get_single_GRF_result_df_para, mp_args_list):
				if df is not None:
					ok_to_keep = (seen_records - total_seen_records)
					total_seen_records = total_seen_records.union(seen_records)
					for seqid in ok_to_keep:
						print(f'{seqid}\t{df[seqid]}', file = out)
				ct += 1
				if ct % 100 == 0:
					print(f'{ct} of {argnum} complete!')
	
	return None



def _get_GRF_result_df_native(genome_name: str, flag_debug: bool) -> Optional[pd.DataFrame]:
	GRF_result_dir_name = f"{genome_name}_GRFmite"
	GRF_result_file_path = os.path.join(GRF_result_dir_name, "candidate.fasta")

	try:
		df_data_dict = [{"id": rec.id, "seq": str(rec.seq), "len": len(rec)}
						for rec in SeqIO.parse(GRF_result_file_path, "fasta")]
	except FileNotFoundError:
		return None

	if not flag_debug:
		subprocess.Popen(["rm", "-rf", GRF_result_dir_name])
	return pd.DataFrame(df_data_dict, columns=["id", "seq", "len"]).astype({"len": int})


def execute(TIRLearner_instance):
	genome_file: str = TIRLearner_instance.genome_file_path
	genome_name: str = TIRLearner_instance.genome_name
	TIR_length: int = TIRLearner_instance.TIR_length
	processors: int = TIRLearner_instance.processors
	para_mode: str = TIRLearner_instance.para_mode
	flag_debug: bool = TIRLearner_instance.flag_debug
	GRF_path: str = TIRLearner_instance.GRF_path

	if NO_PARALLEL in TIRLearner_instance.additional_args:
		return _run_GRF_native(genome_file, genome_name, TIR_length, processors, flag_debug, GRF_path)
	elif para_mode == "gnup":
		raise NotImplementedError()
	else:
		output = os.path.abspath(TIRLearner_instance.processed_de_novo_result_file_name_grf)
		_run_GRF_py_para(genome_file, genome_name, TIR_length, processors, flag_debug, GRF_path,
								TIRLearner_instance.split_fasta_files_path_list, output)
	
		return output
