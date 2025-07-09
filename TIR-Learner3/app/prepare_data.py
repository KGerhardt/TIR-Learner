#!/usr/app/env python3
# -*- coding: utf-8 -*-

from shared import *

def execute(TIRLearner_instance, df_homo: pd.DataFrame = None) -> pd.DataFrame:
	split_size = 500
	chunks_to_process = []
	
	if df_homo is not None:
		remove = set(df_homo["id"].unique())
	
	if "TIRvish" in TIRLearner_instance.keys():
		split = 1
		ct = 0
		next_records = []
		with open(TIRLearner_instance["TIRvish"]) as fh:
			for line in fh:
				next_records.append(line)
				ct += 1
				if ct > 499:
					oname = f'{TIRLearner_instance["TIRvish"]}_chunk_{split}.tsv'
					with open(oname, "w") as out:
						out.write(''.join(next_records))
					next_records = []
					ct = 0
					split += 1
					chunks_to_process.append(os.path.abspath(oname))
						
		if ct > 0:
			oname = f'{TIRLearner_instance["TIRvish"]}_chunk_{split}.tsv'
			with open(oname, "w") as out:
				out.write(''.join(next_records))
			chunks_to_process.append(os.path.abspath(oname))
		
	if "GRF" in TIRLearner_instance.keys():
		split = 1
		ct = 0
		next_records = []
		with open(TIRLearner_instance["GRF"]) as fh:
			for line in fh:
				next_records.append(line)
				ct += 1
				if ct > 499:
					oname = f'{TIRLearner_instance["GRF"]}_chunk_{split}.tsv'
					with open(oname, "w") as out:
						out.write(''.join(next_records))
					next_records = []
					ct = 0
					split += 1
					chunks_to_process.append(os.path.abspath(oname))
						
		if ct > 0:
			oname = f'{TIRLearner_instance["GRF"]}_chunk_{split}.tsv'
			with open(oname, "w") as out:
				out.write(''.join(next_records))
			chunks_to_process.append(os.path.abspath(oname))
	
	return chunks_to_process
