import json
import os

import pyfastx
import re

import numpy as np

seqid_offset_regex = re.compile(r'.+;;(\d+)')
#seqid_seqname_regex = re.compile(r'(.+);;\d+')
seqid_seqname_regex = re.compile(r';;\d+')

#The idea here is good but we need to think of the most useful way to do this all
'''
(1) This is needed by the CNN step, and chunkwise loading via the individual, per-file JSON records is actually a great way to minimize RAM usage
(2) CNN + GRF loads per file should still be modest; passing whatever results to the appropriate CNN model, which searches + creates another JSON of outputs might be best

So the better way to do this all is that the CNN step should 
'''

class json_structure:
	def __init__(self, seqlens, include_label = False):
		self.seqlens = seqlens
		self.offsets = {}
		for s in seqlens:
			offset = 0
			#Check for genomeSplitter structure
			mat = re.match(seqid_offset_regex, s)
			if mat:
				offset = int(mat.groups()[0])
				
			self.offsets[s] = offset
			
		self.include_label = include_label
				
		self.json_record = {}
		
		self.has_records = False
		
	def add_record(self, seqid, start, stop, tsd1, tsd2, tir1, tir2, tir_label = None):
		self.has_records = True
		if seqid not in self.json_record:
			if self.include_label:
				self.json_record[seqid] = {'seq_length':self.seqlens[seqid],
											'chunking_offset':self.offsets[seqid],
											'seq_start_incl_tsd':[],
											'seq_stop_incl_tsd':[],
											'tsd1_size':[],
											'tsd2_size':[],
											'tir1_size':[],
											'tir2_size':[],
											'final_tir_label':[]}
			else:
				self.json_record[seqid] = {'seq_length':self.seqlens[seqid],
											'chunking_offset':self.offsets[seqid],
											'seq_start_incl_tsd':[],
											'seq_stop_incl_tsd':[],
											'tsd1_size':[],
											'tsd2_size':[],
											'tir1_size':[],
											'tir2_size':[]}
		
		self.json_record[seqid]['seq_start_incl_tsd'].append(start)
		self.json_record[seqid]['seq_stop_incl_tsd'].append(stop)
		self.json_record[seqid]['tsd1_size'].append(tsd1)
		self.json_record[seqid]['tsd2_size'].append(tsd2)
		self.json_record[seqid]['tir1_size'].append(tir1)
		self.json_record[seqid]['tir2_size'].append(tir2)
		if self.include_label:
			self.json_record[seqid]['final_tir_label'].append(tir_label)
		
	def sort_records(self):
		#GRF doesn't order starts and stops in its output, so I correct it here because bedtools is faster if the starts are sorted and it's easier to merge if they are, anyway
		for seqid in self.json_record:
			self.json_record[seqid]['seq_start_incl_tsd'] = np.array(self.json_record[seqid]['seq_start_incl_tsd'])
			self.json_record[seqid]['seq_stop_incl_tsd'] = np.array(self.json_record[seqid]['seq_stop_incl_tsd'])
			#The sizes of these should always fit into int16, easier to pass over pickle
			self.json_record[seqid]['tsd1_size'] = np.array(self.json_record[seqid]['tsd1_size'], dtype = np.int16)
			self.json_record[seqid]['tsd2_size'] = np.array(self.json_record[seqid]['tsd2_size'], dtype = np.int16)
			self.json_record[seqid]['tir1_size'] = np.array(self.json_record[seqid]['tir1_size'], dtype = np.int16)
			self.json_record[seqid]['tir2_size'] = np.array(self.json_record[seqid]['tir2_size'], dtype = np.int16)
			
			if self.include_label:
				self.json_record[seqid]['final_tir_label'] = np.array(self.json_record[seqid]['final_tir_label'])
			
			#Order by start (ascending) end (descending), but EXCLUDING the TSD. This causes the longest candidate 
			#at each start locus to be first; at the end of the program this makes resolving overlaps easier
			
			#Numpy lexsort orders in reverse order; the final item (starts) is the primary key, 
			#secondary, tertiary... keys proceed in sort order from right to left
			ordering = np.lexsort((-1 * self.json_record[seqid]['seq_stop_incl_tsd'] - self.json_record[seqid]['tsd2_size'], 
										self.json_record[seqid]['seq_start_incl_tsd'] + self.json_record[seqid]['tsd1_size'],))
			

			self.json_record[seqid]['seq_start_incl_tsd'] = self.json_record[seqid]['seq_start_incl_tsd'][ordering].tolist()
			self.json_record[seqid]['seq_stop_incl_tsd'] = self.json_record[seqid]['seq_stop_incl_tsd'][ordering].tolist()
			self.json_record[seqid]['tsd1_size'] = self.json_record[seqid]['tsd1_size'][ordering].tolist()
			self.json_record[seqid]['tsd2_size'] = self.json_record[seqid]['tsd2_size'][ordering].tolist()
			self.json_record[seqid]['tir1_size'] = self.json_record[seqid]['tir1_size'][ordering].tolist()
			self.json_record[seqid]['tir2_size'] = self.json_record[seqid]['tir2_size'][ordering].tolist()

			if self.include_label:
				self.json_record[seqid]['final_tir_label'] = self.json_record[seqid]['final_tir_label'][ordering].tolist()


#Might need to have an additional version that decomposes TIRvish + GRF records into their parts for CNN work; this is OK I think
class json_loader:
	def __init__(self, working_dir = None):
		self.json_data = None
		self.workloads = None
		self.sequence_names = None
		
		self.work_dir = working_dir
		#Very unlikely to see this pattern unless it's made by genomeSplitter
		self.path_regex = re.compile(r'.+(split_genome.+_chunk_.+_offset_.+\.fasta)')
		
	def load_json(self, json_file, get_names = True):
		with open(json_file) as inf:
			self.json_data = json.load(inf)
			
		self.workloads = []
		
		if get_names:
			self.sequence_names = []
		else:
			self.sequence_names = None
			
		for source_file in self.json_data:
			new_source = source_file
			if self.work_dir is not None:
				if self.work_dir not in source_file:
					base = re.match(self.path_regex, source_file)
					if base:
						base = base.groups()[0]
						new_source = os.path.join(self.work_dir, base)
						#print(f'From existing GRF or TIRvish JSON, source updated from {source_file} to {new_source}')
					else:
						print(f'{source_file} not found, this will provbably fail.')
			
			this_workload = {}
			for seqid in self.json_data[source_file]:
				if len(self.json_data[source_file]) > 0:
					if get_names:
						this_workload[seqid] = {'starts':[],'ends':[],'names':[]}
					else:
						this_workload[seqid] = {'starts':[],'ends':[],'names':None}

					for s, e, tsd1, tsd2, tir1, tir2 in	zip(
								self.json_data[source_file][seqid]['seq_start_incl_tsd'], self.json_data[source_file][seqid]['seq_stop_incl_tsd'],
								self.json_data[source_file][seqid]['tsd1_size'], self.json_data[source_file][seqid]['tsd2_size'],
								self.json_data[source_file][seqid]['tir1_size'], self.json_data[source_file][seqid]['tir2_size']
								):
						
						this_workload[seqid]['starts'].append(s-1)
						this_workload[seqid]['ends'].append(e)
						
						if get_names:
							name = f'{seqid}:start={s}:stop={e}:tsd1_sz={tsd1}:tsd2_sz={tsd2}:tir1_sz={tir1}:tir2_sz={tir2}'
							this_workload[seqid]['names'].append(name)
			
			if len(this_workload) > 0:
				self.workloads.append((new_source, this_workload,))
	
	def load_json_for_cnn(self, json_file):
		with open(json_file) as inf:
			self.json_data = json.load(inf)
			
		self.workloads = []
		for source_file in self.json_data:
			new_source = source_file
			if self.work_dir is not None:
				if self.work_dir not in source_file:
					base = re.match(self.path_regex, source_file)
					if base:
						base = base.groups()[0]
						#Update to new path
						new_source = os.path.join(self.work_dir, base)
						#print(f'From existing GRF or TIRvish JSON, source updated from {source_file} to {new_source}')
					else:
						print(f'{source_file} not found, this will provbably fail.')
						
			self.workloads.append((new_source, self.json_data[source_file],))
	
#class for converting JSON records to bed files, calling bedtools, working with those files, etc.
class bed_worker:
	def __init__(self, json, has_names = True, working_directory = '.'):
		self.wd = working_directory
		
		self.json_chunk = json
		self.has_names = has_names
		
		self.source = None
		
		self.my_ref_genome = None
		self.my_seqlens = None
		self.my_loaded_sequences = None
		
		self.homolog_filter_indices = None
		
		#No longer needed
		self.bed_file_name = None
		self.bed_file_contents = None
		
		self.seqid_metadata = None
		self.reverse_json_record = None
	
	def load_refgen(self):
		self.source = self.json_chunk[0]
		self.my_ref_genome = {}
		self.my_seqlens = {}
		
		for name, seq in pyfastx.Fasta(self.source, build_index = False):
			self.my_ref_genome[name] = seq
			self.my_seqlens[name] = len(seq)
		
	#This is much faster with direct in-mem slicing via pyfastx + indices
	def convert_json_to_sequences(self, index_names = False):
		self.my_loaded_sequences = {}
		
		if index_names:
			self.homolog_filter_indices = {}
		
		for seqid in self.json_chunk[1]:
			if index_names:
				homolog_index = 0
				self.homolog_filter_indices[seqid] = {}
				
			if self.has_names:
				for s, e, n in zip(self.json_chunk[1][seqid]['starts'], self.json_chunk[1][seqid]['ends'], self.json_chunk[1][seqid]['names']):
					self.my_loaded_sequences[n] = self.my_ref_genome[seqid][s:e]
					if index_names:
						self.homolog_filter_indices[seqid][n] = homolog_index
						homolog_index += 1
					
			else:
				for s, e in zip(self.json_chunk[1][seqid]['starts'], self.json_chunk[1][seqid]['ends']):
					n = f'{seqid}:start={s}:stop={e}'
					self.my_loaded_sequences[n] = self.my_ref_genome[seqid][s:e]
					if index_names:
						self.homolog_filter_indices[seqid][n] = homolog_index
						homolog_index += 1
	
	#This is much faster with direct in-mem slicing via pyfastx + indices;
	#The minimum size here is to prohibit checking sequences that TIR-Learner considers too short 
	def convert_json_to_sequences_for_cnn(self, minimum_seq_size = 50):
		self.my_loaded_sequences = []
		
		seq_index = 0
		self.seqid_indices = {}
		self.reverse_json_record = {}
		self.seqid_metadata = {}
		
		for seqid in self.json_chunk[1]:
			self.seqid_metadata[seqid] = (self.json_chunk[1][seqid]['seq_length'], self.json_chunk[1][seqid]['chunking_offset'],)
			for s, e, tsd1, tsd2, tir1, tir2 in	zip(
							self.json_chunk[1][seqid]['seq_start_incl_tsd'], self.json_chunk[1][seqid]['seq_stop_incl_tsd'],
							self.json_chunk[1][seqid]['tsd1_size'], self.json_chunk[1][seqid]['tsd2_size'],
							self.json_chunk[1][seqid]['tir1_size'], self.json_chunk[1][seqid]['tir2_size']
							):
				
				#The full string with TSDs
				this_subsequence = self.my_ref_genome[seqid][(s-1):e]
				#Size of the sequence without TSDs
				this_seqlen = len(this_subsequence) - (tsd1 + tsd2)
				if this_seqlen >= minimum_seq_size:
					
					#Already has TSDs subtracted, so only needs TIRs removed now
					middle_section_size = this_seqlen - (tir1 + tir2)
					
					start_loc = 0
					tsd1_string = this_subsequence[start_loc:(start_loc + tsd1)]
					start_loc += tsd1
					tir1_string = this_subsequence[start_loc:(start_loc + tir1)]
					start_loc += tir1
					
					#Need to get the middle string part
					if middle_section_size > 0:
						middle_string = this_subsequence[start_loc:(start_loc + middle_section_size)]
						start_loc += middle_section_size
					else:
						middle_string = ''
					
					tir2_string = this_subsequence[start_loc:(start_loc + tir2)]
					start_loc += tir2
					tsd2_string = this_subsequence[start_loc:(start_loc + tsd2)]
					#start_loc += tsd1
					
					digested_entry = (tsd1_string, tir1_string, middle_string, tir2_string, tsd2_string,)
					
					reverse_record = (seqid, s, e, tsd1, tsd2, tir1, tir2, )
					#print(tsd1_string + tir1_string + middle_string + tir2_string + tsd2_string == this_subsequence)
					
					self.my_loaded_sequences.append(digested_entry)
					self.reverse_json_record[seq_index] = reverse_record
					
					#Only increment these if a sequence is included / keep the seq index effectively 
					#0:(num_records -1) so that numpy where as a proxy index works later
					seq_index += 1
	
	def convert_json_to_sequences_for_BLAST(self, tsd_pcts, tir_pcts, minimum_seq_size = 0,):
		self.my_loaded_sequences = []
		
		tir_types = []
		clean_tsd_pcts = []
		clean_tir_pcts = []
		
		seq_index = 0
		self.seqid_indices = {}
		self.reverse_json_record = {}
		self.seqid_metadata = {}
		
		for seqid in self.json_chunk[1]:
			self.seqid_metadata[seqid] = (self.json_chunk[1][seqid]['seq_length'], self.json_chunk[1][seqid]['chunking_offset'],)
			for s, e, tsd1, tsd2, tir1, tir2, tir_type, tsd_pct, tir_pct in zip(
							self.json_chunk[1][seqid]['seq_start_incl_tsd'], self.json_chunk[1][seqid]['seq_stop_incl_tsd'],
							self.json_chunk[1][seqid]['tsd1_size'], self.json_chunk[1][seqid]['tsd2_size'],
							self.json_chunk[1][seqid]['tir1_size'], self.json_chunk[1][seqid]['tir2_size'],
							self.json_chunk[1][seqid]['final_tir_label'], tsd_pcts[seqid], tir_pcts[seqid]
							):
				
				#The full string with TSDs
				this_subsequence = self.my_ref_genome[seqid][(s-1):e]
				#Size of the sequence without TSDs
				this_seqlen = len(this_subsequence) - (tsd1 + tsd2)
				if this_seqlen >= minimum_seq_size:
					#Already has TSDs subtracted, so only needs TIRs removed now
					middle_section_size = this_seqlen - (tir1 + tir2)
					
					start_loc = 0
					tsd1_string = this_subsequence[start_loc:(start_loc + tsd1)]
					start_loc += tsd1
					tir1_string = this_subsequence[start_loc:(start_loc + tir1)]
					start_loc += tir1
					
					#Need to get the middle string part
					if middle_section_size > 0:
						middle_string = this_subsequence[start_loc:(start_loc + middle_section_size)]
						start_loc += middle_section_size
					else:
						middle_string = ''
					
					tir2_string = this_subsequence[start_loc:(start_loc + tir2)]
					start_loc += tir2
					tsd2_string = this_subsequence[start_loc:(start_loc + tsd2)]
					#start_loc += tsd1
					
					digested_entry = (tsd1_string, tir1_string, middle_string, tir2_string, tsd2_string,)
					
					reverse_record = (seqid, s, e, tsd1, tsd2, tir1, tir2, )
					
					self.my_loaded_sequences.append(digested_entry)
					self.reverse_json_record[seq_index] = reverse_record
					
					tir_types.append(tir_type)
					clean_tsd_pcts.append(tsd_pct)
					clean_tir_pcts.append(tir_pct)
					
					#Only increment these if a sequence is included / keep the seq index effectively 
					#0:(num_records -1) so that numpy where as a proxy index works later
					seq_index += 1
					
					
		return tir_types, clean_tsd_pcts, clean_tir_pcts
	
	def cnn_filter_json(self, passing_indices, tir_types, tsd_percents, tir_percents, module = 'Module4', cnn_scores = None):				
		new_json = {}
		finalized_sequence = {}
		finalized_gff3 = {}
		
		'''
		#GFF3 format line looks like:
		Chr1	Module4	DTM	39311	39451	.	.	.	TIR:GGGGCAATTG_CAAATGCCCC_90.0_TSD:AAAAAATAA_AAAAAATAA_100.0-+-141
		
		#FASTA format line looks like
		>TIR-Learner_Chr1_39311_39451_DTM_TIR:GGGGCAATTG_CAAATGCCCC_90.0_TSD:AAAAAATAA_AAAAAATAA_100.0-+-141
		[GGGGCAATTG]CAAATTTACCACTGCTTTGAATCGTTATTGCAAGGCTACCACTAGAAAATTGCAAAAATGCCACTAGAACTGTCA...
		
		#Fasta starts with TIR, does NOT incl TSD

		'''
		
		if cnn_scores is None:
			cnn_scores = [None] * len(passing_indices)
		
		for index, tir_type, tsd_pct, tir_pct, cnn in zip(passing_indices, tir_types, tsd_percents, tir_percents, cnn_scores):
			#Abbreviated names for tir and tsd sizes
			seqid, s, e, ts1, ts2, ti1, ti2 = self.reverse_json_record[index]
			
			stripped_seqname = re.sub(seqid_seqname_regex, '', seqid)
			if seqid not in new_json:
				finalized_sequence[seqid] = []
				finalized_gff3[seqid] = []
				
				new_json[seqid] = {'seq_length':self.seqid_metadata[seqid][0],
											'chunking_offset':self.seqid_metadata[seqid][1],
											'seq_start_incl_tsd':[],
											'seq_stop_incl_tsd':[],
											'tsd1_size':[],
											'tsd2_size':[],
											'tsd_percent':[],
											'tir1_size':[],
											'tir2_size':[],
											'tir_percent':[],
											'final_tir_label':[]}
				
				if cnn is not None:
					new_json[seqid]['cnn_scores'] = []
			
			new_json[seqid]['seq_start_incl_tsd'].append(s)
			new_json[seqid]['seq_stop_incl_tsd'].append(e)
			new_json[seqid]['tsd1_size'].append(ts1)
			new_json[seqid]['tsd2_size'].append(ts2)
			new_json[seqid]['tsd_percent'].append(tsd_pct)
			new_json[seqid]['tir1_size'].append(ti1)
			new_json[seqid]['tir2_size'].append(ti2)
			new_json[seqid]['tir_percent'].append(tir_pct)
			new_json[seqid]['final_tir_label'].append(tir_type)
			if cnn is not None:
				new_json[seqid]['cnn_scores'].append(cnn)
			
			off = new_json[seqid]['chunking_offset']
			
			tir_only_seq = ''.join(self.my_loaded_sequences[index][1:4])
			
			ele_size = len(tir_only_seq)
			tir1 = self.my_loaded_sequences[index][1]
			tir2 = self.my_loaded_sequences[index][3]
			tsd1 = self.my_loaded_sequences[index][0]
			tsd2 = self.my_loaded_sequences[index][4]
			
			if cnn is None:
				suffix = f'TIR:{tir1}_{tir2}_{tir_pct}_TSD:{tsd1}_{tsd2}_{tsd_pct}-+-{ele_size}'
			else:
				cnn_label = '\t'.join([str(s) for s in cnn])
				suffix = f'TIR:{tir1}_{tir2}_{tir_pct}_TSD:{tsd1}_{tsd2}_{tsd_pct}-+-{ele_size}\t{cnn_label}'
			
			#The final records need the offset-corrected values
			one_idx_start_no_tsd = s + ts1 + off
			one_idx_stop_no_tsd = e - ts2 + off
			
			gff_record = '\t'.join([stripped_seqname,
							module,
							tir_type,
							str(one_idx_start_no_tsd),
							str(one_idx_stop_no_tsd),
							'.',
							'.',
							'.',
							suffix])
							
			finalized_gff3[seqid].append(gff_record)
			
			#name = f'{seqid}:start={s}:stop={e}:tsd1_sz={tsd1}:tsd2_sz={tsd2}:tir1_sz={tir1}:tir2_sz={tir2}'
			final_seqname = f'>TIR-Learner_{stripped_seqname}_{one_idx_start_no_tsd}_{one_idx_stop_no_tsd}_{tir_type}_{suffix}'

			#Encode these are one record for sanity later
			finalized_sequence[seqid].append(f'{final_seqname}\n{tir_only_seq}')

		#Resolve overlaps
		keep_indices = {}
		for seqid in new_json:
			'''
			
			The desired goal here is to remove all overlaps except a full, internal nest.
			
			A, B = TIR candidates
			m = TIR start or end locus
			- = internal base
			y = internal element (after TIR 1 end, before TIR 2 start) start or end locus
			= = TIR base
			
			This is the ONLY kind of overlap we permit, where the entire 
			TIR structure of B fits within only the internal portion of A
			y1=======m1------------A------------m2=======y2
			             y1===m1---B---m2===y2
						 
			All other overlaps are resolved by removing the shorter candidate
			
			'''
			
			#Earlier in the program the json_structure.sort_records() function sorts by start ascending and stop descending;
			#this means that the first appearance of a record left by this point is always the longest element at that start index
			
			
			
			#Collect starts + ends without TSDs
			starts = np.array(new_json[seqid]['seq_start_incl_tsd']) + np.array(new_json[seqid]['tsd1_size'])
			ends   = np.array(new_json[seqid]['seq_stop_incl_tsd']) -  np.array(new_json[seqid]['tsd2_size'])
			element_only_starts = starts + np.array(new_json[seqid]['tir1_size'])
			element_only_ends   = ends   - np.array(new_json[seqid]['tir2_size'])
			
			keep_indices[seqid] = np.zeros(starts.shape[0], dtype = np.bool_)
			
			#Remove repeat start indices; keep the first index of each, only, 
			#because that's the longest candidate at that start due to sort order
			#These kinds of overlaps are removed by definition as tir-tir overlaps; retain longest only
			mask = np.ones(starts.shape[0], dtype=bool)
			mask[1:] = starts[1:] != starts[:-1]
			starts = starts[mask]
			ends = ends[mask]
			element_only_starts = element_only_starts[mask]
			element_only_ends = element_only_ends[mask]
			
			#These are actually 1 less than the element size but the comparison is all that matters
			element_sizes = ends - starts
			
			#This mask corresponds to the indices of the pre-filtered json arrays; 
			#it also corresponds to entries in finalized_gff3 and finalized_sequence arrays per seqid for easy filt
			mask_loci = np.where(mask)
			
			current_size = starts.shape[0]
			
			#If an end is greater than the next start, that end and its corresponding start are an overlap pair
			has_tir_overlap = np.where(ends[:-1] >= starts[1:])[0]
			
			'''
			Okay, what the hell do I need to do?
			
			same start elements have already been cleaned by the above...
			
			Maybe we need to iterate remove only bad overlaps until the size of the arrays doesn't change, 
			skipping over nests
			'''
			
			#kicked = 0
			#rounds = 0
			#nests = 0
			#We may need to repeatedly iterate;
			#loop ends here if there are no remaining overlaps at all
			while has_tir_overlap.shape[0] > 0:
				#rounds += 1
				#nests = 0
				
				remove_indices = set()
				for i in has_tir_overlap:
					#If we have already removed an element that overlaps with the next element, it does not bear checking
					if i not in remove_indices:
						#Find all elements which this TIR candidate overlaps; adjust their indices to the overall list
						my_overlaps = np.where(np.logical_and(starts[i+1:] > starts[i], starts[i+1:] <= ends[i]))[0] + i + 1
						for j in my_overlaps:
							#if an overlap is nested; skip it
							if element_only_starts[i] < starts[j] and element_only_ends[i] > ends[j]:
								pass
								#print('Nesting element found:')
								#print(f'{element_only_starts[i]} to {element_only_ends[i]}')
								#print('NESTS') 
								#print(f'{starts[j]} to {ends[j]}')
								#nests += 1
							else:
								#Original favored long sequences, but tir-learner's long sequences tend to be worse; try reversing it?
								#if element_sizes[i] > element_sizes[j]:
								if element_sizes[i] > element_sizes[j]:
									#print(f'Element i {i} kicks element j {j}')
									#We are removing element j
									remove_indices.add(i)
								else:
									#print(f'Element j {j} kicks element i {i} and another loop might happen')
									#The current element is going to be removed, which means information
									#about what else it overlaps is no longer of any value and we break the loop;
									#other possible overlaps will be handled by more iterations
									remove_indices.add(j)
									break
									
				#There are still elements left to remove / cleaning is not yet done
				if len(remove_indices) > 0:
					#kicked += len(remove_indices)
					remove_indices = np.array(list(remove_indices))
					
					#Keep track of the original loci for filtering other categories later on...
					mask_loci = np.delete(mask_loci, remove_indices)
					
					#Remove the loci that are no longer needed
					starts = np.delete(starts, remove_indices)
					ends = np.delete(ends, remove_indices)
					element_only_starts = np.delete(element_only_starts, remove_indices)
					element_only_ends = np.delete(element_only_ends, remove_indices)
					
					#If there are any nested overlaps, this will have length > 0 even if there is nothing remaining to remove
					has_tir_overlap = np.where(ends[:-1] >= starts[1:])[0]
				
				#Exit the loop if there are only nested elements
				else:
					break
		
			#Store the loci of final sequences to be kept
			keep_indices[seqid][mask_loci] = True
		
		return new_json, finalized_gff3, finalized_sequence, keep_indices
			
	def fake_fasta(self):
		self.my_loaded_sequences = ''.join(f'>{k}\n{v}\n' for k, v in self.my_loaded_sequences.items())
	
#5200 is the default TIR size limit (5k) + 200 (extension size default, not really used anymore)
def dereplicate_json(json_data, overlap_size = 5200):
	chunk_id_regex = re.compile(r'long_chunk_(.+)_offset_(\d+).fasta')
	
	cleaned_json = {}
	
	#Purge repeats as we load this data to prevent redundant CNN effort, outputs
	groupings = {}
	for k in json_data.keys():
		#Huh. So this doesn't guarantee that every index for this grouping will be present. 
		#It will only be the ones in the json, which could include skips
		if 'long_chunk_' in k:
			if len(json_data[k]) > 0:
				cleaned_json[k] = {}
				for v in json_data[k]:
					cleaned_json[k][v] = {'seq_length':json_data[k][v]['seq_length'],
										'chunking_offset':json_data[k][v]['chunking_offset'],
										'seq_start_incl_tsd':[],
										'seq_stop_incl_tsd':[],
										'tsd1_size':[],
										'tsd2_size':[],
										'tir1_size':[],
										'tir2_size':[]}
				
				#I think this os.path.basename and the TIRvish stuff are fighting in some genomes
				mat = re.match(chunk_id_regex, os.path.basename(k)).groups()
				seqid = mat[0]
			
			else:
				print(f'JSON key {k} had no entries and will be skipped.')
			
			off = int(mat[1])
			if seqid not in groupings:
				groupings[seqid] = []
				
			groupings[seqid].append((k, off,))
		else:
			if len(json_data[k]) > 0:
				#Manually copy short groups to the new JSON, otherwise with dicts it's done by reference
				cleaned_json[k] = {}
				for v in json_data[k]:
					cleaned_json[k][v] = {'seq_length':json_data[k][v]['seq_length'],
										'chunking_offset':json_data[k][v]['chunking_offset'],
										'seq_start_incl_tsd':json_data[k][v]['seq_start_incl_tsd'],
										'seq_stop_incl_tsd':json_data[k][v]['seq_stop_incl_tsd'],
										'tsd1_size':json_data[k][v]['tsd1_size'],
										'tsd2_size':json_data[k][v]['tsd2_size'],
										'tir1_size':json_data[k][v]['tir1_size'],
										'tir2_size':json_data[k][v]['tir2_size']}
				json_data[k] = None
				
		
	#The logic of this section is to look at long chunk groups and sort them by offset from 0 to whatever
	#The idea is to look at the next section and see what items overlap in that section with the tail of the current section
	#Purge the current section's tail and add it to the cleaned json
	#Then move to the next chunk
	#The final chunk of each long sequence is always added intact
	
	for seqid in groupings:
		groupings[seqid] = sorted(groupings[seqid], key=lambda x: x[1])
		
		for i in range(0, len(groupings[seqid]) - 1):
			src, off = groupings[seqid][i]
			
			next_src, next_off = groupings[seqid][i+1]
						
			indices_to_remove = []
			
			#cutoff = next_off - overlap_size
			
			access_name = f'{seqid};;{off}'
			next_access_name = f'{seqid};;{next_off}'
			
			#We could probably do some math to figure out where exactly the sequnces 
			#overlap because that should be all possible repeats
			#But it's just not much work to brute force it here
			
			these_starts = None
			these_ends = None
			
			next_starts = None
			next_ends = None
			
			if src in json_data:
				if access_name in json_data[src]:
					these_starts = np.array(json_data[src][access_name]['seq_start_incl_tsd']) + off
					these_ends   = np.array(json_data[src][access_name]['seq_stop_incl_tsd']) + off
					if next_src in json_data:
						if next_access_name in json_data[next_src]:								
							next_starts = np.array(json_data[next_src][next_access_name]['seq_start_incl_tsd']) + next_off
							next_ends = np.array(json_data[next_src][next_access_name]['seq_stop_incl_tsd']) + next_off
							
							#Find shared start indices in the genome, get their indicess
							intersect_values, left_starts, right_starts = np.intersect1d(these_starts, next_starts, return_indices = True)
							#Often, there will be no intersections
							if left_starts.shape[0] > 0:
								#Subset to only the overlaps
								next_starts, next_ends = next_starts[right_starts], next_ends[right_starts]
								#Double loop is dangerous but it's very small arrays in basically all cases
								for i in left_starts:
									this_start = these_starts[i]
									this_end = these_ends[i]
									for j in next_ends[next_starts == this_start]:
										if this_end == j:
											indices_to_remove.append(int(i))
											
								indices_to_remove = np.unique(np.array(indices_to_remove, dtype = np.int64))
					
					valid_indices = np.arange(len(these_starts))
					#If there are skip indices, remove those
					if len(indices_to_remove) > 0:
						valid_indices = np.delete(valid_indices, indices_to_remove)
					
					#Copy the data over from the current long chunk
					for i in valid_indices:
						cleaned_json[src][access_name]['seq_start_incl_tsd'].append(json_data[src][access_name]['seq_start_incl_tsd'][i])
						cleaned_json[src][access_name]['seq_stop_incl_tsd'].append(json_data[src][access_name]['seq_stop_incl_tsd'][i])
						cleaned_json[src][access_name]['tsd1_size'].append(json_data[src][access_name]['tsd1_size'][i])
						cleaned_json[src][access_name]['tsd2_size'].append(json_data[src][access_name]['tsd2_size'][i])
						cleaned_json[src][access_name]['tir1_size'].append(json_data[src][access_name]['tir1_size'][i])
						cleaned_json[src][access_name]['tir2_size'].append(json_data[src][access_name]['tir2_size'][i])
					
					#Free space
					json_data[src] = None
				else:
					print(f'Access name {access_name} not in json_data[{src}]. This should not happen. Please report the error.')
			else:
				print(f'Src {src} not in json data. This should not happen. Please report the error.')				
					
		#Final group			
		src, off = groupings[seqid][len(groupings[seqid]) - 1]
		access_name = f'{seqid};;{off}'
		if src in json_data:
			if access_name in json_data[src]:
				for i in range(0, len(json_data[src][access_name]['seq_start_incl_tsd'])):
					cleaned_json[src][access_name]['seq_start_incl_tsd'].append(json_data[src][access_name]['seq_start_incl_tsd'][i])
					cleaned_json[src][access_name]['seq_stop_incl_tsd'].append(json_data[src][access_name]['seq_stop_incl_tsd'][i])
					cleaned_json[src][access_name]['tsd1_size'].append(json_data[src][access_name]['tsd1_size'][i])
					cleaned_json[src][access_name]['tsd2_size'].append(json_data[src][access_name]['tsd2_size'][i])
					cleaned_json[src][access_name]['tir1_size'].append(json_data[src][access_name]['tir1_size'][i])
					cleaned_json[src][access_name]['tir2_size'].append(json_data[src][access_name]['tir2_size'][i])
			else:
				print(f'Access name {access_name} not in json_data[{src}]. This should not happen. Please report the error.')
		else:
			print(f'Src {src} not in json data. This should not happen. Please report the error.')
	
	#Final pass, remove anything empty
	final_keys = list(cleaned_json.keys())
	for k in final_keys:
		if len(cleaned_json[k]) == 0:
			remove = cleaned_json.pop(k)
		else:
			for jk in list(cleaned_json[k].keys()):
				if len(cleaned_json[k][jk]['tsd1_size']) == 0:
					remove = cleaned_json[k].pop(jk)
	
	return cleaned_json

'''
#Testing json dereplication code
jj = 'TIRVish_json.txt'
with open(jj) as inf:
	dat = json.load(inf)

new_j = dereplicate_json(dat)
'''

'''	
#This works as a template function for others to call and work with easily	
def parabed(load):
	bl = bed_worker(load, has_names = get_names, working_directory = work_dir)
	bl.load_refgen()
	bl.convert_json_to_sequences()
	
	#This conversion basically works as an in-memory string rep of the fasta file to pass to something like BLAST
	
	#Probably need a... pandas? equivalent for CNN calls. IDK, see what it needs
	payload = ''.join(f'>{k}\n{v}\n' for k, v in bl.my_loaded_sequences.items())
	with open(os.path.join(work_dir, f'{os.path.basename(bl.source)}.membed.fasta'), 'w') as out:
		out.write(payload)
	
	return bl.source
	

mn = json_loader()

import sys
f = sys.argv[1]

global work_dir
work_dir = 'TIR_Learner_working_directory/current_results'

global get_names

get_ext = False
get_names = True

mn.load_json(f, retrieve_extended = get_ext, get_names = get_names)
#mn.load_json(f, retrieve_extended = True)

with multiprocessing.Pool(20) as pool:
	for r in pool.imap_unordered(parabed, mn.workloads):
		pass

#mn.collect()
'''