#!/usr/app/env python3
# -*- coding: utf-8 -*-

from shared import *


def _get_sequence_fragment(x: pd.Series, feature_size: int = 200) -> str:
	seq = x["seq"]
	len_seq = len(seq)
	if len_seq >= feature_size * 2:
		return seq[0:feature_size] + seq[-feature_size:]
	s1 = seq[0:int(len_seq / 2)]
	s2 = seq[int(len_seq / 2):]
	n1 = "N" * (feature_size - len(s1))
	n2 = "N" * (feature_size - len(s2))
	s1 += n1
	s2 = n2 + s2
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

def run_by_chunks(chunk):
	df = pd.read_csv(chunk, sep = "\t", dtype = {'id':str, 'seq':str}, names = ("id", "seq"))

	#print("  Step 1/7: Getting sequence fragment for prediction")
	df["seq_frag"] = df.swifter.progress_bar(False).apply(_get_sequence_fragment, axis=1)
	#just process in-flow rather than loading the data again later
	df = df.drop(columns="seq")

	df = _feature_encoding(df, False)
	df = _ml_predict(df, CNN_MODEL_DIR_ABS_PATH)
	if df is not None:
		df = _post_processing(df, False)
	
	os.remove(chunk)
	
	return df

	
#partially from get_fasta_sequence
def _get_start_end(fasta_len_dict : dict, df_in: pd.DataFrame, flag_verbose: bool = True) -> pd.DataFrame:
	df = df_in.copy()
	df["start"] = df.swifter.progress_bar(False).apply(lambda x: min(x["sstart"], x["send"]), axis=1)
	df["end"] = df.swifter.progress_bar(False).apply(lambda x: max(x["sstart"], x["send"]), axis=1)
	df["sstart"] = df["start"]
	df["send"] = df["end"]

	df["start"] = df.swifter.progress_bar(False).apply(lambda x: max(x["sstart"] - TE_TOLERANCE, 1), axis=1)
	# start = start - length if start > length else 1
	df["end"] = df.swifter.progress_bar(False).apply(lambda x: x["send"] + TE_TOLERANCE, axis=1)

	# Ensure "end" not exceeding seq's length
	df["end"] = df.swifter.progress_bar(False).apply(lambda x: min(x["end"], fasta_len_dict[x["seqid"]]), axis=1)

	#df = df.drop_duplicates(["start", "end", "seqid", "TIR_type"], keep="first", ignore_index=True)
	return df

def execute(TIRLearner_instance) -> pd.DataFrame:
	#num_procs = TIRLearner_instance.processors
	#chunks = TIRLearner_instance["base"]
	#We actually want to pass a filepath here and load chunks
	#df = TIRLearner_instance["base"].copy()

	num_chunks = min([TIRLearner_instance.processors, len(TIRLearner_instance["base"])])
	#gl = TIRLearner_instance.genome_lengths
	
	#dflist = []
	
	ct = 0
	print(f"{len(TIRLearner_instance["base"])} partitions to process with CNN")
	pool = mp.Pool(num_chunks)
	with open(TIRLearner_instance.processed_de_novo_result_file_name_cnn, "w", newline='') as out:
		for result in pool.imap_unordered(run_by_chunks, TIRLearner_instance["base"]):
			ct += 1
			if result is not None:
				result = _get_start_end(TIRLearner_instance.genome_lengths, result)
				#for col in result.columns:
				#	print(f"Column: {col}, Data Type: {result[col].dtype}")
				
				#dflist.append(result)
				result.to_csv(out, sep = '\t', header = False, index = False)
			print(f'{ct} of {num_chunks} complete!')
		
	#df = pd.concat(dflist).reset_index(drop=True)

	return TIRLearner_instance.processed_de_novo_result_file_name_cnn