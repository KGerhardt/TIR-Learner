#!/usr/app/env python3
# -*- coding: utf-8 -*-

from shared import *

from get_fasta_sequence import get_fasta_pieces_SeqIO_post
from process_de_novo_result import TA_repeats_check


def _combine_all(filepaths) -> Optional[pd.DataFrame]:
	
	'''
	if len(df_list) > 1:
		try:
			df = pd.concat(df_list, ignore_index=True)
		except ValueError:
			df = None
	else:
		df = df_list[0].copy()
	'''
	
	
	'''
	pd.read_csv(TIRLearner_instance["base"], sep = "\t", dtype = {'TIR_type':str, 'id':str, 'seqid':str, 'sstart':int, 'send':int, 'start':int, 'end':int, 'seq':str}, 
										names = ('TIR_type', 'id', 'seqid', 'sstart', 'send', 'start', 'end', 'seq'),
										chunksize = 100_000):
	'''
	df_list = []
	for f in filepaths:
		df = pd.read_csv(filepaths, sep = '\t', 
						dtype = {"seqid":str, "source":str, "type":str, "sstart":int, "send":int, "TIR1":str, "TIR2":str, "TIR_percent":float, "TSD1":str, "TSD2":str, "TSD_percent":float, "len":int},
						names = ("seqid", "source", "type", "sstart", "send", "TIR1", "TIR2", "TIR_percent", "TSD1", "TSD2", "TSD_percent", "len"))
		df_list.append(df)
	
	if len(df_list) > 1:
		try:
			df = pd.concat(df_list, ignore_index=True)
		except ValueError:
			df = None
	else:
		df = df_list[0].copy()
	
	if df is None or df.shape[0] == 0:
		return None

	df = df.sort_values(["seqid", "sstart", "send", "source", "type"], ignore_index=True)
	df = df.drop_duplicates(["seqid", "sstart", "send"], keep="first", ignore_index=True)
	df = df.drop_duplicates(["seqid", "sstart"], keep="first", ignore_index=True)
	df = df.drop_duplicates(["seqid", "send"], keep="first", ignore_index=True)

	df["TIR_pair_str"] = df["TIR1"] + df["TIR2"]
	df = TA_repeats_check(df, "TIR_pair_str")

	if df.shape[0] == 0:
		return None

	return df.drop(columns="TIR_pair_str")


def _format_df_in_gff3_format(df_in: pd.DataFrame, flag_verbose: bool = True) -> pd.DataFrame:
	df = df_in.copy()
	df["attributes"] = df.swifter.progress_bar(flag_verbose).apply(
		lambda x: (f"TIR:{x['TIR1']}_{x['TIR2']}_{x['TIR_percent']}_"
				   f"TSD:{x['TSD1']}_{x['TSD2']}_{x['TSD_percent']}{SPLITER}{x['len']}"), axis=1)
	df = df.loc[:, ["seqid", "source", "type", "sstart", "send", "attributes"]]
	df.insert(5, "phase", ".")
	df.insert(5, "strand", ".")
	df.insert(5, "score", ".")
	return df


# =============================== Remove the Shorter One Among Two Overlapped Sequences ================================

# noinspection PyUnusedLocal
def check_element_overlap(x1: int, y1: int, x2: int, y2: int) -> bool:
	"""Checking sequences overlap only among element part of sequences

	Precondition: x1 < x2, y1 != y2

	  - True: Overlap, x2 <= y1 < y2

		x1------A------y1
			   x2------B------y2

	  - False: Nesting, y2 < y1

		x1--------A---------y1
			  x2----B----y2

	  - False: No correlation, y1 < x2

		x1----A----y1
					   x2------B------y2
	"""
	if x2 <= y1 < y2:
		return True
	return False

# noinspection PyUnusedLocal
def check_element_TIR_overlap(x1: int, y1: int, x2: int, y2: int, m1: int, n1: int, m2: int, n2: int) -> bool:
	"""Checking sequences overlap among element part and TIR part of sequences

	Precondition: x1 < x2, y1 < x2 or y1 > y2

	  - True: Overlap - A right element with B left TIR

		m1===x1------A------y1===n1
					   m2========x2------B------y2========n2

	  - True: Overlap - B left element with A right TIR

		m1========x1------A------y1========n1
								 m2===x2------B------y2===n2

	  - True: Overlap - A right element with B right TIR

		m1===x1----------------A----------------y1===n1
				  m2========x2------B------y2========n2

	  - True: Overlap - A left element with B left TIR

		m1===x1----------------A----------------y1===n1
		m2========x2------B------y2========n2

	  - False: No correlation

		m1===x1----A----y1===n1
								 m2=====x2------B------y2=====n2
	"""
	if y1 < x2:
		if y1 > m2 or x2 < n1:
			return True
	elif y1 > y2:
		if y1 < n2 or x1 > m2:
			return True
	return False


# def remove_overlap(df_in):
#	 df = df_in.sort_values(by=["sstart", "send"], ignore_index=True)
#	 df["TIR_len"] = df.swifter.progress_bar(True).apply(lambda x: len(x["TIR"][0]), axis=1)
#	 df["tstart"] = df.loc[:, "sstart"] - df.loc[:, "TIR_len"]
#	 df["tend"] = df.loc[:, "send"] + df.loc[:, "TIR_len"]
#	 idx_sstart = df.columns.get_loc("sstart")
#	 idx_send = df.columns.get_loc("send")
#	 idx_tstart = df.columns.get_loc("tstart")
#	 idx_tend = df.columns.get_loc("tend")
#	 ptr = 0
#	 while ptr+1 < df.shape[0]:
#		 if (check_element_overlap(df.iloc[ptr, idx_sstart], df.iloc[ptr, idx_send],
#								   df.iloc[ptr+1, idx_sstart], df.iloc[ptr+1, idx_send]) or
#				 check_element_TIR_overlap(df.iloc[ptr, idx_sstart], df.iloc[ptr, idx_send],
#										   df.iloc[ptr+1, idx_sstart], df.iloc[ptr+1, idx_send],
#										   df.iloc[ptr, idx_tstart], df.iloc[ptr, idx_tend],
#										   df.iloc[ptr+1, idx_tstart], df.iloc[ptr+1, idx_tend])):
#			 df = df.drop(labels=df.iloc[[ptr, ptr+1], df.columns.get_loc("len")].idxmin()).reset_index(drop=True)
#		 ptr += 1
#	 df = df.drop(columns=["TIR_len", "tstart", "tend"])
#	 return df


def _remove_overlap(df_in: pd.DataFrame, flag_verbose: bool = True) -> pd.DataFrame:
	"""
	TODO documentation needed
	:param df_in:
	:param flag_verbose:
	:return:
	"""
	df = df_in.sort_values(by=["sstart", "send"], ignore_index=True)
	seqid = df.loc[0, "seqid"]
	df["TIR_len"] = df["TIR1"].str.len()
	df["tstart"] = df.loc[:, "sstart"] - df.loc[:, "TIR_len"]
	df["tend"] = df.loc[:, "send"] + df.loc[:, "TIR_len"]
	idx_len = df.columns.get_loc("len")
	idx_sstart = df.columns.get_loc("sstart")
	idx_send = df.columns.get_loc("send")
	idx_tstart = df.columns.get_loc("tstart")
	idx_tend = df.columns.get_loc("tend")

	dropped_index_list = []
	ptr1 = 0
	ptr2 = 1
	while True:
		while ptr1 in dropped_index_list:
			ptr1 += 1
		while ptr1 >= ptr2 or ptr2 in dropped_index_list:
			ptr2 += 1
		if not ptr1 < ptr2 < df.shape[0]:
			break
		if (check_element_overlap(df.iloc[ptr1, idx_sstart], df.iloc[ptr1, idx_send],
								  df.iloc[ptr2, idx_sstart], df.iloc[ptr2, idx_send]) or
				check_element_TIR_overlap(df.iloc[ptr1, idx_sstart], df.iloc[ptr1, idx_send],
										  df.iloc[ptr2, idx_sstart], df.iloc[ptr2, idx_send],
										  df.iloc[ptr1, idx_tstart], df.iloc[ptr1, idx_tend],
										  df.iloc[ptr2, idx_tstart], df.iloc[ptr2, idx_tend])):
			dropped_index_list.append(int(df.iloc[[ptr1, ptr2], idx_len].idxmin()))
		ptr1 += 1
		ptr2 += 1
	df = df.drop(dropped_index_list)
	df = df.drop(columns=["TIR_len", "tstart", "tend"])
	if dropped_index_list:
		if flag_verbose:
			print(f"Sequence(s) removed from genome {seqid} (index): {dropped_index_list}")
		else:
			print(f"{len(dropped_index_list)} sequence(s) of genome {seqid} are removed")
	return df

# ======================================================================================================================


def _get_final_fasta_file(df_in: pd.DataFrame, genome_file: str, genome_name: str, processors: int, file_path: str,
						  flag_verbose: bool = True):
	df = df_in.copy()
	df["name"] = df.swifter.progress_bar(flag_verbose).apply(
		lambda x: f">{genome_name}_{x['seqid']}_{x['sstart']}_{x['send']}_{x['type']}_{x['attributes']}", axis=1)
	df.rename(columns={"sstart": "start", "send": "end"}, inplace=True)
	
	df = get_fasta_pieces_SeqIO_post(genome_file, df, file_path, processors, flag_verbose)
	
	#df = df.loc[:, ["name", "seq"]]
	#df.to_csv(file_path, index=False, header=False, sep="\n")


def execute(TIRLearner_instance, raw_result_df_list: List[pd.DataFrame]):
	genome_file: str = TIRLearner_instance.genome_file_path
	genome_name: str = TIRLearner_instance.genome_name
	result_output_dir_path: str = os.path.join(TIRLearner_instance.output_dir_path, RESULT_OUTPUT_DIR_NAME)
	processors: int = TIRLearner_instance.processors
	flag_verbose: bool = TIRLearner_instance.flag_verbose

	terminal_print("Post Processing")

	'''
	df = df.sort_values(["seqid", "sstart", "send", "source", "type"], ignore_index=True)
	df = df.drop_duplicates(["seqid", "sstart", "send"], keep="first", ignore_index=True)
	df = df.drop_duplicates(["seqid", "sstart"], keep="first", ignore_index=True)
	df = df.drop_duplicates(["seqid", "send"], keep="first", ignore_index=True)

	df["TIR_pair_str"] = df["TIR1"] + df["TIR2"]
	df = TA_repeats_check(df, "TIR_pair_str")
	'''

	print("  Step 1/6: Combining all results")
	df_combined = _combine_all(TIRLearner_instance['m4'])
			
	if df_combined is None:
		print("[WARN] No TIR found. Post-processing will be terminated and no result will be produced.")
		return

	os.makedirs(result_output_dir_path, exist_ok=True)

	print("  Step 2/6: Preparing gff3 attributes for all sequences")
	df_gff3 = _format_df_in_gff3_format(df_combined, flag_verbose)
	df_gff3.to_csv(os.path.join(result_output_dir_path, f"{genome_name}_FinalAnn.gff3"), index=False, header=False,
				   sep="\t")

	print("  Step 3/6: Generating raw fasta file")
	_get_final_fasta_file(df_gff3, genome_file, genome_name, processors,
						  os.path.join(result_output_dir_path, f"{genome_name}_FinalAnn.fa"), flag_verbose)
	del df_gff3
	gc.collect()

	df_combined_groupby_seqid = df_combined.groupby("seqid")
	df_combined_seqid_list = [df_combined_groupby_seqid.get_group(df) for df in df_combined_groupby_seqid.groups]
	df_combined_mp = [(df, flag_verbose) for df in df_combined_seqid_list]
	del df_combined, df_combined_groupby_seqid, df_combined_seqid_list

	print("  Step 4/6: Removing the shorter one among two overlapped sequences")
	if not flag_verbose:
		print("	", end="")
	with mp.Pool(processors) as pool:
		df_filtered_list = pool.starmap(_remove_overlap, df_combined_mp)
	if not flag_verbose:
		print('\n', end="")
	df_filtered = pd.concat(df_filtered_list, ignore_index=True)
	del df_filtered_list

	print("  Step 5/6: Preparing gff3 attributes for all sequences")
	df_gff3_filtered = _format_df_in_gff3_format(df_filtered, flag_verbose)
	df_gff3_filtered.to_csv(os.path.join(result_output_dir_path, f"{genome_name}_FinalAnn_filter.gff3"), index=False,
							header=False, sep="\t")

	print("  Step 6/6: Generating final fasta file")
	_get_final_fasta_file(df_gff3_filtered, genome_file, genome_name, processors,
						  os.path.join(result_output_dir_path, f"{genome_name}_FinalAnn_filter.fa"), flag_verbose)
