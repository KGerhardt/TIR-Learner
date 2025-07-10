#!/usr/app/env python3
# -*- coding: utf-8 -*-

from shared import *

from get_fasta_sequence import get_fasta_pieces_SeqIO_tirvish


# Use noqa to suppress "Shadowing built-ins" -------------------↴
def create_sequence_record(seq: str, id: str) -> SeqRecord:  # noqa
	"""
	Create a BioPython SeqRecord object from a sequence string and ID.

	This function takes a biological sequence string and an identifier,
	and returns a SeqRecord object from the BioPython library.

	:param seq: The biological sequence as a string
	:type seq: str
	:param id: The unique identifier for the sequence
	:type id: str
	:return: A SeqRecord object containing the sequence and its metadata
	:rtype: SeqRecord

	Example::

		>>> record = create_sequence_record("ATGC", "seq1")
		>>> print(record.id)
		seq1
		>>> print(record.seq)
		ATGC
	"""
	return SeqRecord(Seq(seq), id=id, description="")


def split_sequence_evenly(seq_record: SeqRecord, split_seq_len: int, overlap_seq_len: int) -> List[SeqRecord]:
	"""
	Evenly splits a sequence into several segmented sequences, each segment has a maximum length of
	the length threshold split_seq_len.

	Parameters:
	- seq_record: A SeqRecord object representing the sequence.
	- split_seq_len: The length threshold for splitting, i.e. maximum allowed sequence length before splitting.
	- overlap_seq_len: An overlapping segment that covers the boundary of the preceding segment and
	the succeeding segment with the length of overlap_seq_len will be created. overlap_seq_len // 2 represents
	the length of the two parts that come from the preceding segment and the succeeding segment.
	Specifying a value <=1 will result in no overlapping segment being produced.

	Returns:
	A list of SeqRecord objects representing the original or split sequences.
	"""
	records = []
	seq_len = len(seq_record)
	overlap_seg_half_len = overlap_seq_len // 2
	total_parts = seq_len // split_seq_len + (1 if seq_len % split_seq_len else 0)
	for i in range(total_parts):
		segment_start = i * split_seq_len
		segment_end = min(segment_start + split_seq_len, seq_len)
		segment_seq = str(seq_record.seq[segment_start:segment_end])
		part_id = MP_SPLIT_SEQ_ID_FORMAT_STR.format(seq_record.id, i + 1, total_parts)
		# part_id = f"{seq_record.id}_split_{i + 1}of{total_parts}"
		records.append(create_sequence_record(segment_seq, part_id))
		# records.append(SeqRecord(segment_seq, id=part_id, description=""))
		if overlap_seg_half_len > 0 and segment_end != seq_len:
			overlap_seg = str(seq_record.seq[segment_end - overlap_seg_half_len:
										 min(segment_end + overlap_seg_half_len, seq_len)])
			overlap_id = MP_SPLIT_SEQ_ID_FORMAT_STR.format(seq_record.id, f"{i + 1}.5", total_parts)
			# overlap_id = f"{seq_record.id}_split_{i + 1}.5of{total_parts}"
			records.append(create_sequence_record(overlap_seg, overlap_id))
			# records.append(SeqRecord(overlap_seg, id=overlap_id, description=""))
	return records


def _save_fasta_file(seq_record: SeqRecord) -> str:
	fasta_file_name = f"{seq_record.id}.fasta"
	working_dir_name = f"{fasta_file_name}_{SPLIT_FASTA_TAG}"
	os.makedirs(working_dir_name, exist_ok=True)
	fasta_file_path = os.path.join(working_dir_name, fasta_file_name)
	SeqIO.write(seq_record, fasta_file_path, "fasta")
	return fasta_file_path


def process_fasta(genome_file: str, split_seq_len: int, overlap_seq_len: int) -> List[str]:
	"""
	Write each sequence in the FASTA file into separate FASTA files and further split the sequence into segments if
	when needed based on a length threshold split_seq_len.

	Parameters:
	- file_name: Path to the FASTA file.
	- split_seq_len: Length threshold for splitting sequence.
	"""
	split_fasta_files_path_list = []
	for seq_record in SeqIO.parse(genome_file, "fasta"):
		if len(seq_record) >= MP_SPLIT_SEQ_LEN:
			segments = split_sequence_evenly(seq_record, split_seq_len, overlap_seq_len)
			for segment in segments:
				split_fasta_files_path_list.append(_save_fasta_file(segment))
		else:
			split_fasta_files_path_list.append(_save_fasta_file(seq_record))
	return split_fasta_files_path_list


def retrieve_split_sequence_offset(segment_position: str, split_seq_len: int, overlap_seq_len: int) -> int:
	if split_seq_len == 0:
		raise ValueError("When split happens, split_seq_len cannot be zero.")
	try:
		# Normal segment, segment_position is like: <x>of<y>, where x and y are integers
		segment_index = int(segment_position.split("of")[0])
		offset = (segment_index - 1) * split_seq_len
	except ValueError:
		# Overlap segment, segment_position is like: <x>.5of<y>, where x and y are integers
		segment_index = int(segment_position.split("of")[0][:-2])  # remove the .5
		offset = segment_index * split_seq_len - overlap_seq_len // 2
	return offset


def _TIRvish(genome_file: str, genome_name: str, TIR_length: int, gt_path: str) -> str:
	gt_bin_path = os.path.join(gt_path, "gt")
	gt_index_file_name = genome_name + SPLITER + "gt_index"
	subprocess.Popen(
		[gt_bin_path, "suffixerator", "-db", genome_file, "-indexname", gt_index_file_name,
		 "-tis", "-suf", "-lcp", "-des", "-ssp", "-sds", "-dna", "-mirrored"]).wait()

	TIRvish_result_gff3_file_name = f"{genome_name}{SPLITER}TIRvish.gff3"
	gt_tirvish = (f"\"{gt_bin_path}\" tirvish -index {gt_index_file_name} -seed 20 -mintirlen 10 -maxtirlen 1000 "
				  f"-mintirdist 10 -maxtirdist {str(TIR_length)} -similar 80 -mintsd 2 -maxtsd 11 "
				  f"-vic 13 -seqids \"yes\" > {TIRvish_result_gff3_file_name}")
	subprocess.Popen(gt_tirvish, shell=True).wait()
	subprocess.Popen(["find", ".", "-name", f"{gt_index_file_name}*", "-delete"])
	return TIRvish_result_gff3_file_name


def _TIRvish_mp(genome_file_path: str, genome_name: str, TIR_length: int, gt_path: str) -> str:
	TIRvish_working_dir = os.path.dirname(genome_file_path)
	genome_file_name = os.path.basename(genome_file_path)
	os.chdir(TIRvish_working_dir)
	TIRvish_result_gff3_file_name = _TIRvish(genome_file_name, genome_name, TIR_length, gt_path)
	os.chdir("../")
	return os.path.join(TIRvish_working_dir, TIRvish_result_gff3_file_name)


def _get_TIRvish_result_df(TIRvish_result_gff3_file_path: str, flag_debug: bool,
						   split_seq_len: int = 0, overlap_seq_len: int = 0) -> pd.DataFrame:
	df_data_dict = {"seqid": [], "start": [], "end": [], "TIR1_start": [], "TIR1_end": [], "TIR2_start": [],
					"TIR2_end": [], "id": []}
	df_type = {"start": int, "end": int, "TIR1_start": int, "TIR1_end": int, "TIR2_start": int, "TIR2_end": int}

	if os.path.exists(TIRvish_result_gff3_file_path) and os.path.getsize(TIRvish_result_gff3_file_path) != 0:
		with open(TIRvish_result_gff3_file_path, 'r') as f:
			while (line := f.readline()) != "":
				while line.startswith('#'):
					line = f.readline()
				TSD1 = list(map(int, f.readline().split('\t')[3:5]))
				record = f.readline().split('\t')
				TIR1 = list(map(int, f.readline().split('\t')[3:5]))
				TIR2 = list(map(int, f.readline().split('\t')[3:5]))
				TSD2 = list(map(int, f.readline().split('\t')[3:5]))
				f.readline()  # Jump the line "###"

				seqid = record[0]
				start = int(record[3])
				end = int(record[4])

				offset = 0
				# Offset for segment of sequence from splitting
				if MP_SPLIT_SEQ_TAG in seqid:
					# seqid = <original_seqid><MP_SPLIT_SEQ_TAG><segment_position>
					segment_position = seqid.split(MP_SPLIT_SEQ_TAG)[1]
					seqid = seqid.split(MP_SPLIT_SEQ_TAG)[0]  # reassign original_seqid to seqid
					offset = retrieve_split_sequence_offset(segment_position, split_seq_len, overlap_seq_len)

				# -1 at start coordinates are internal adjustments for get_fasta_sequence of TIR-Learner
				# Actual start coordinates are without the -1 adjustments, i.e. start + offset
				df_data_dict["seqid"].append(seqid)
				df_data_dict["start"].append(start - 1 + offset)
				df_data_dict["end"].append(end + offset)
				df_data_dict["TIR1_start"].append(TIR1[0] - 1 + offset)
				df_data_dict["TIR1_end"].append(TIR1[1] + offset)
				df_data_dict["TIR2_start"].append(TIR2[0] - 1 + offset)
				df_data_dict["TIR2_end"].append(TIR2[1] + offset)
				df_data_dict["id"].append(
					f">{seqid}:{start + offset}:{end + offset}:"
					f"tsd1_{TSD1[0] + offset}_{TSD1[1] + offset}_tsd2_{TSD2[0] + offset}_{TSD2[1] + offset}_"
					f"tir1_{TIR1[0] + offset}_{TIR1[1] + offset}_tir2_{TIR2[0] + offset}_{TIR2[1] + offset}")

	if not flag_debug:
		subprocess.Popen(["unlink", TIRvish_result_gff3_file_path])
		
	return pd.DataFrame(df_data_dict).astype(df_type)


def _run_TIRvish_native(genome_file: str, genome_name: str, TIR_length: int,
						flag_debug: bool, gt_path: str) -> pd.DataFrame:
	print("  Step 1/2: Executing TIRvish in native mode")
	TIRvish_result_gff3_file_name = _TIRvish(genome_file, genome_name, TIR_length, gt_path)
	print("  Step 2/2: Getting TIRvish result")
	return _get_TIRvish_result_df(TIRvish_result_gff3_file_name, flag_debug)


def _run_TIRvish_py_para(genome_file: str, genome_name: str, TIR_length: int,
						 processors: int, flag_debug: bool, gt_path: str,
						 fasta_files_path_list: List[str]) -> pd.DataFrame:
						 
						 
	os.makedirs(f"{SPLIT_FASTA_TAG}_mp", exist_ok=True)
	os.chdir(f"./{SPLIT_FASTA_TAG}_mp")

	print("  Step 1/3: Processing FASTA files")
	fasta_files_path_list.extend(process_fasta(genome_file, MP_SPLIT_SEQ_LEN, MP_OVERLAP_SEQ_LEN))

	print("  Step 2/3: Executing TIRvish with python multiprocessing")
	mp_args_list = [(file_path, genome_name, TIR_length, gt_path) for file_path in fasta_files_path_list]
	with mp.Pool(processors) as pool:
		TIRvish_result_gff3_file_path_list = pool.starmap(_TIRvish_mp, mp_args_list)

	print("  Step 3/3: Getting TIRvish result")
	mp_args_list = [(file_path, flag_debug, MP_SPLIT_SEQ_LEN, MP_OVERLAP_SEQ_LEN) for file_path in
					TIRvish_result_gff3_file_path_list]
					
	with mp.Pool(processors) as pool:
		TIRvish_result_df_list = pool.starmap(_get_TIRvish_result_df, mp_args_list)

	os.chdir("../")
	return pd.concat(TIRvish_result_df_list).reset_index(drop=True)


def execute(TIRLearner_instance) -> pd.DataFrame:
	genome_file = TIRLearner_instance.genome_file_path
	genome_name = TIRLearner_instance.genome_name
	TIR_length = TIRLearner_instance.TIR_length
	processors = TIRLearner_instance.processors
	para_mode = TIRLearner_instance.para_mode
	# TODO add GNU Parallel support
	flag_verbose = TIRLearner_instance.flag_verbose
	flag_debug = TIRLearner_instance.flag_debug
	gt_path = TIRLearner_instance.gt_path
	additional_args = TIRLearner_instance.additional_args
	fasta_files_path_list = TIRLearner_instance.split_fasta_files_path_list

	if NO_PARALLEL in additional_args:
		df = _run_TIRvish_native(genome_file, genome_name, TIR_length, flag_debug, gt_path)
	elif para_mode == "gnup":
		raise NotImplementedError()
	else:
		#This code is still not low-mem'd, but that should be OK
		df = _run_TIRvish_py_para(genome_file, genome_name, TIR_length, processors, flag_debug, gt_path,
								  fasta_files_path_list)
								  
	#Here is where cleaning goes:
	df = df[df["end"] - df["start"] + 1 >= 50].copy() #initial length clean
	#Add TIR loc values
	df["TIR1_start"] = df["TIR1_start"] - df["start"]
	df.loc[df["TIR1_start"] < 0, "TIR1_start"] = 0
	df["TIR1_end"] = df["TIR1_end"] - df["start"]

	df["TIR2_start"] = df["TIR2_start"] - df["start"]
	df.loc[df["TIR2_start"] < 0, "TIR2_start"] = 0
	df["TIR2_end"] = df["TIR2_end"] - df["start"]
	
	output_file = TIRLearner_instance["TIRvish"]
	
	#Cleaning has been pushed to the sequence retrieval component so that sequences which are not needed are not kept
	df = get_fasta_pieces_SeqIO_tirvish(genome_file, df, output_file, processors, flag_verbose)

	#Final cleanup to only the columns process_de_novo_results would have produced
	#df = df.dropna(ignore_index=True).loc[:, ["id", "seq"]].copy()

	return None
