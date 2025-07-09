#!/usr/app/env python3
# -*- coding: utf-8 -*-

from shared import *

#Pushed cleaning to this function so that unneeded seqs are not loaded by tirvish
def clean_record(seqid, t1s, t1e, t2s, t2e, sequence):
	is_ok = False
	tir1 = sequence[t1s:t1e]
	tir2 = sequence[t2s:t2e]

	seqlen = len(sequence)
	
	max_N_pct = 0.2
	max_ta_pct = 0.7
	
	cts = Counter(sequence)
	#TA and N checks on sequence
	if (cts['N'] / seqlen) < max_N_pct and ((cts['A']+cts['T']) / seqlen) < max_ta_pct:
		cts = Counter(tir1)
		tirlen = len(tir1)
		#TA and N checks on tir1
		if (cts['N'] / tirlen) < max_N_pct and ((cts['A']+cts['T']) / tirlen) < max_ta_pct:
			cts = Counter(tir2)
			tirlen = len(tir2)
			#TA and N checks on tir2
			if (cts['N'] / tirlen) < max_N_pct and ((cts['A']+cts['T']) / tirlen) < max_ta_pct:
				#TIRvish does not include a TSD check like GRF, so this is the last check
				is_ok = True
	
	if is_ok:
		return sequence
	else:
		return pd.NA

def get_fasta_pieces_single_seqid_SeqIO_tirvish(args) -> Optional[pd.DataFrame]:
	df_in, genome_SeqRecord, seqid, flag_verbose = args

	df = df_in[df_in["seqid"] == seqid].copy()
	if df.shape[0] == 0:
		return None

	# [IMPORTANT] Indexing Conversion
	# 1-based indexing seq slice [start, end] <=> 0-based indexing Python slice [start-1, end) = [start-1, end-1]
	df["seq"] = df.swifter.progress_bar(flag_verbose).apply(
		lambda x: clean_record(x["seqid"], x["TIR1_start"], x["TIR1_end"], x["TIR2_start"], x["TIR2_end"], str(genome_SeqRecord.seq[x["start"] - 1: x["end"]])), axis=1)
	
	#Clean the record
	df = df.dropna(ignore_index=True).loc[:, ["id", "seq"]].copy()
	
	return df
	
def get_fasta_pieces_SeqIO_tirvish(genome_file: str, df_in: pd.DataFrame, output_file:str,
						   processors: int, flag_verbose: bool = True) -> Optional[pd.DataFrame]:
	# if not df_in cannot be used as it will cause "ValueError: The truth value of a DataFrame is ambiguous."
	if df_in is None or df_in.shape[0] == 0:
		return None

	df = df_in.copy()
	
	#I see, it's indexing the fasta file and then passing the whole dataframe to each processor to load sequences for that seqid, then joining them at the end
	genome_SeqIO_index = SeqIO.index(genome_file, "fasta")
	mp_args_list = [(df, genome_SeqIO_index[seqid], seqid, flag_verbose) for seqid in genome_SeqIO_index]
	
	with open(output_file, "w", newline = '') as out:
		with mp.Pool(processors) as pool:
			for result in pool.imap_unordered(get_fasta_pieces_single_seqid_SeqIO_tirvish, mp_args_list):
				result.to_csv(out, sep = '\t', header = False, index = False)
		
	return None


#Original functions, unchanged for use later on
def get_fasta_pieces_single_seqid_SeqIO(df_in: pd.DataFrame, genome_SeqRecord: SeqIO.SeqRecord, seqid: str,
										flag_verbose: bool = True) -> Optional[pd.DataFrame]:
	df = df_in[df_in["seqid"] == seqid].copy()
	if df.shape[0] == 0:
		return None

	# [IMPORTANT] Indexing Conversion
	# 1-based indexing seq slice [start, end] <=> 0-based indexing Python slice [start-1, end) = [start-1, end-1]
	df["seq"] = df.swifter.progress_bar(flag_verbose).apply(
		lambda x: str(genome_SeqRecord.seq[x["start"] - 1: x["end"]]), axis=1)
		
	df = df.dropna()
	
	return df

def get_fasta_pieces_SeqIO(genome_file: str, df_in: pd.DataFrame,
						   processors: int, flag_verbose: bool = True) -> Optional[pd.DataFrame]:
	# if not df_in cannot be used as it will cause "ValueError: The truth value of a DataFrame is ambiguous."
	if df_in is None or df_in.shape[0] == 0:
		return None

	df = df_in.copy()
	
	#I see, it's indexing the fasta file and then passing the whole dataframe to each processor to load sequences for that seqid, then joining them at the end
	genome_SeqIO_index = SeqIO.index(genome_file, "fasta")
	mp_args_list = [(df, genome_SeqIO_index[seqid], seqid, flag_verbose) for seqid in genome_SeqIO_index]
	
	with mp.Pool(processors) as pool:
		df_with_seq_list = pool.starmap(get_fasta_pieces_single_seqid_SeqIO, mp_args_list)
		
	return pd.concat(df_with_seq_list).sort_index()


def _get_start_end(genome_file: str, df_in: pd.DataFrame, tolerance: int, flag_verbose: bool = True) -> pd.DataFrame:
	df = df_in.copy()
	df["start"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: min(x["sstart"], x["send"]), axis=1)
	df["end"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: max(x["sstart"], x["send"]), axis=1)
	df["sstart"] = df["start"]
	df["send"] = df["end"]

	df["start"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: max(x["sstart"] - tolerance, 1), axis=1)
	# start = start - length if start > length else 1
	df["end"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: x["send"] + tolerance, axis=1)

	# Ensure "end" not exceeding seq's length
	fasta_len_dict = {rec.id: len(rec.seq) for rec in SeqIO.parse(genome_file, "fasta")}
	df["end"] = df.swifter.progress_bar(flag_verbose).apply(lambda x: min(x["end"], fasta_len_dict[x["seqid"]]), axis=1)

	df = df.drop_duplicates(["start", "end", "seqid", "TIR_type"], keep="first", ignore_index=True)
	return df


def execute(TIRLearner_instance) -> pd.DataFrame:
	df = _get_start_end(TIRLearner_instance.genome_file_path, TIRLearner_instance["base"], TE_TOLERANCE,
						TIRLearner_instance.flag_verbose)
	return get_fasta_pieces_SeqIO(TIRLearner_instance.genome_file_path, df, TIRLearner_instance.processors,
								  TIRLearner_instance.flag_verbose)
