import sys
import os
import multiprocessing
import subprocess
import json

import shutil
import re

import pyfastx

from .genomeSplitter import genomeSplitter
from .new_seq_reader import json_loader, bed_worker, json_structure
from .new_tir_tsd import tsd_tir_checker
from .output_compressor import compress

import numpy as np

program_root = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

reflib_regex = re.compile(r'(.+)##(.+)##(.+)')

def make_one_db(infile):
	output_database = os.path.join(workdir, 'module1', f'{os.path.basename(infile)}.blast.db')
	output_database_check = os.path.join(workdir, 'module1', f'{os.path.basename(infile)}.blast.db.complete')
	
	if not os.path.exists(output_database_check):
		comm = f'makeblastdb -in {infile} -out {output_database} -parse_seqids -dbtype nucl'
		comm = comm.split()
		
		subprocess.run(comm, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
		with open(output_database_check, 'w') as out:
			pass
			
	return output_database


def blast_ref_vs_genome(args):
	query, database = args[0], args[1]
	original_sequence = os.path.basename(database).replace('.blast.db', '')
	original_sequence = os.path.join(workdir, 'module1', original_sequence)
	
	blast_format = '6 qseqid sseqid length pident gaps mismatch qstart qend sstart send evalue qcovhsp'
	blast_command = f'blastn -max_hsps 5 -perc_identity 80 -qcov_hsp_perc 100 -query {query} -db {database} -outfmt'
	blast_command = blast_command.split()
	blast_command.append(blast_format)
	
	#print(f'Querying {os.path.basename(query)} vs {os.path.basename(database)}')
	proc = subprocess.run(blast_command, capture_output = True, text=True)
	
	seen_values = {}
	cleaned_proc = {}
	
	#Clean blast output to select one best hit per genomic region
	for line in proc.stdout.splitlines():
		segs = line.strip().split('\t')
		genome_seqid = segs[1]
		
		reference_library_seqid = segs[0]
		
		#Get tir type info about the sequence
		tag = re.match(reflib_regex, reference_library_seqid).groups()
		tir_type = tag[1]
		
		#Skip non-tirs; this really shouldn't ever happen with the ref homolog search
		if tir_type != 'NonTIR':
			#Start + stop, 1-indexed because BLAST
			q1, q2 = int(segs[8]), int(segs[9])
			if q2 > q1:
				start, end = q1, q2
			else:
				start, end = q2, q1
				
			#Skip short sequences
			sub_length =  end - start + 1
			if sub_length >= 50:
				#Hit quality score
				evalue = float(segs[10])
				
				#If the genome seqid is new, add a whole new record, which is a dict of starts: {end : info}
				if genome_seqid not in seen_values:
					seen_values[genome_seqid]  = {start:{end:evalue}}
					cleaned_proc[genome_seqid] = {start:{end:tir_type}}
				else:
					#If it's a new start, it must be a new start + end pair, add it
					if start not in seen_values[genome_seqid]:
						seen_values[genome_seqid][start]  = {end:evalue}
						cleaned_proc[genome_seqid][start] = {end:tir_type}
					else:
						#If it's NOT a new end, compare e values and retain the lowest e-value / the single best match
						if end in seen_values[genome_seqid][start]:
							if evalue < seen_values[genome_seqid][start][end]:
								seen_values[genome_seqid][start][end]  = evalue
								cleaned_proc[genome_seqid][start][end] = tir_type
						else:
							#New end just gets added
							seen_values[genome_seqid][start][end]  = evalue
							cleaned_proc[genome_seqid][start][end] = tir_type
		
	#Clean up an unneeded part
	seen_values = None
	
	#translator for revcomp
	rc = str.maketrans('ATCG', 'TAGC')
	
	#Load the relevant genome chunk
	my_seqs = {}
	seqlens = {}
	for name, seq in pyfastx.Fasta(original_sequence, build_index = False):
		my_seqs[name] = seq		
		seqlens[name] = len(seq)
	
	#Prepare TIR+TSD checker, final output data
	tt = tsd_tir_checker()
	tsd_check_size = 10
	
	result_json = json_structure(seqlens, include_label = True)
	
	tsd_pcts = {}
	tir_pcts = {}
	#process cleaned BLAST output for TIR and TSD presence/absence; most other checks have already been done
	for genome_seqid in cleaned_proc:
		for start in cleaned_proc[genome_seqid]:
			for end in cleaned_proc[genome_seqid][start]:
				tir_type = cleaned_proc[genome_seqid][start][end]
		
				#This is the in-genome homology TIR sequence we're looking at - it should be from TIR-TIR but probably not include TSD
				substring = my_seqs[genome_seqid][start-1:end]
				
				#Because we don't know how big the TIR will be, we extract the longest valid de-novo TIR
				midpt = int(sub_length / 2)
				#Break the sequence in half
				front = substring[0:midpt]
				#Select back half + reverse + complement
				back = substring[-midpt:][::-1].translate(rc)
				
				#Align the halves and check for a >80% similar aligned invert
				has_tir, left_tir_size, right_tir_size, ref_starts_at, query_starts_at, pct = tt.wfa_align(front, back, min_size = 10, min_similarity = 0.8)
					
				#If such a region is found, proceed
				if has_tir:
					#Check if the TIR starts or ends with the correct sequence for the recovered type
					conserved_ok = tt.check_tir_conservation(tir_type, front, back)
					
					#If the TIR start is correct, proceed
					if conserved_ok:
						
						#Extend the genome slice 10 bp left and right
						tsd_left_start = max([0, start - tsd_check_size-1])
						left_tsd_string = my_seqs[genome_seqid][tsd_left_start:(tsd_left_start+tsd_check_size)]
						right_tsd_string = my_seqs[genome_seqid][(end):(end+tsd_check_size-1)]
						
						#Find the longest correctly positioned TSD in that region at or above 80% seqid
						has_tsd, left_tsd_size, right_tsd_size, tsd_percent = tt.check_tsd(left_tsd_string, 
																			right_tsd_string, 
																			tir_type = tir_type,
																			min_similarity = 0.8)
						
						#If such a TSD exists, add it as a record
						if has_tsd:
							#Stops with TSD included
							full_element_start = start - left_tsd_size
							full_element_stop = end + right_tsd_size
							
							if genome_seqid not in tsd_pcts:
								tsd_pcts[genome_seqid] = []
								tir_pcts[genome_seqid] = []
							
							tsd_pcts[genome_seqid].append(tsd_percent)
							tir_pcts[genome_seqid].append(pct)
							
							#This should include tsd and TIR percentages...
							result_json.add_record(seqid = genome_seqid, 
												start = full_element_start, 
												stop = full_element_stop, 
												tsd1 = left_tsd_size, 
												tsd2 = right_tsd_size,
												tir1 = left_tir_size, 
												tir2 = right_tir_size, 
												tir_label = tir_type)
	
	#Sort records for use
	result_json.sort_records()
	
	#Spoof loading the genome since we already did that once
	#Python copy by ref for dicts is actually useful, here
	bl = bed_worker((original_sequence, result_json.json_record), has_names = True)
	bl.source = original_sequence
	bl.my_ref_genome = my_seqs
	bl.my_seqlens = seqlens
	
	tir_types, tsd_pcts, tir_pcts = bl.convert_json_to_sequences_for_BLAST(tsd_pcts, tir_pcts, minimum_seq_size = 0)
	#for this code, passing indices is just 1:num_seqs
	clean_json, final_gff3, final_fasta, keep_indices = bl.cnn_filter_json(passing_indices = list(
																							range(0, len(bl.my_loaded_sequences))
																							),
																			tir_types = tir_types,
																			tsd_percents = tsd_pcts,
																			tir_percents = tir_pcts,
																			module = 'Module1')
																			
	
	return original_sequence, clean_json, final_gff3, final_fasta, keep_indices
	
#This is the code for executing a de-novo sequences vs. reference libraries search.
#
def json_blast(args):
	workload, target_database = args
	bl = bed_worker(workload, has_names = True, working_directory = '.')
	
	bl.load_refgen()
	bl.convert_json_to_sequences(index_names = True)
	bl.fake_fasta()
	
	blast_format = '6 qseqid sseqid length pident gaps mismatch qstart qend sstart send evalue qcovhsp'
	blast_command = f'blastn -max_hsps 5 -perc_identity 80 -qcov_hsp_perc 80 -db {target_database} -outfmt'
	
	blast_command = blast_command.split()
	blast_command.append(blast_format)
		
	proc = subprocess.run(blast_command, input = bl.my_loaded_sequences, capture_output = True, text = True)
	
	#All we need to do here is filter the input JSON and write a "no_homologs" output, 
	#which really just means return "no_homolog" indices to main
	
	#Record all the appearances of every unique GRF or TIRvish record that passes 80% id, 80 qcovhsp filters
	sequences_with_homolog_hits = {}
	for line in proc.stdout.splitlines():
		grf_or_tirvish_name = line.strip().split('\t')[0]
		seqid = grf_or_tirvish_name.split(':')[0]
		if seqid not in sequences_with_homolog_hits:
			sequences_with_homolog_hits[seqid] = set([grf_or_tirvish_name])
		else:
			sequences_with_homolog_hits[seqid].add(grf_or_tirvish_name)
	
	#prepare output repository
	non_homolog_records = {}
	for seqid in sequences_with_homolog_hits:
		#Collect all the original sequence names
		all_seqs = set(list(bl.homolog_filter_indices[seqid].keys()))
		
		#Use set difference to find sequences that have no blast hits
		non_homologs = all_seqs - sequences_with_homolog_hits[seqid]
		
		#If there are any, record the positional indices of those sequences within the JSON under that seqid header
		if len(non_homologs) > 0:
			non_homolog_records[seqid] = []
			for key in non_homologs:
				non_homolog_records[seqid].append(bl.homolog_filter_indices[seqid][key])
			
			#And sort those indices for easier retrieval, since key sorting will almost certainly not preserve order
			non_homolog_records[seqid].sort()
			
		#print(bl.source, seqid)
		#print(non_homolog_records[seqid])
		
	return bl.source, non_homolog_records

class blaster:
	def __init__(self, reference_genome, species, working_dir, threads = 1):
		self.rg = reference_genome
		self.species = species
		
		self.species_libraries = None
		
		self.wd = working_dir
		self.blast_dir = os.path.join(self.wd, 'module1')
		if not os.path.exists(self.blast_dir):
			os.makedirs(self.blast_dir, exist_ok = True)
		
		self.threads = threads
		
		self.gs = genomeSplitter(genome_file = self.rg,
								output_directory = self.blast_dir,
								chunk_size = 0,
								#Use a safe overlap size, Ensure that extension is possible later
								overlap_size = 0,
								minimum_seq_size = 0,
								procs = self.threads,
								smart = False,
								post_index = False,
								verbose = True,
								overwrite = False,
								do_bedtools_prep = True)
		
		self.reference_genome_chunks = None
		
		self.ref_blast_databases = None
		
		self.library_blast_dbs = None
		
		self.get_species_ref()
		
		global workdir
		workdir = self.wd
		global ext_size
		ext_size = 200
		
		self.ref_homolog_file = None
		self.tirvish_homologs = None
		self.grf_homologs = None
		
		
	def get_species_ref(self):
		#global species_ref_dict
		#species_ref_dict = {}
	
		reference_libraries = os.path.join(program_root, 'JointRefLib')
		self.species_libraries = []
		for f in os.listdir(reference_libraries):
			if self.species in f:
				if f.endswith('_TEs.fasta'):
					this_library = os.path.join(reference_libraries, f)
					self.species_libraries.append(this_library)
					#for record in pyfastx.Fasta(this_library, build_index = True):
					#	species_ref_dict[record.name] = record.description
			
		#Sort by size to push longer searches earlier
		self.species_libraries.sort(key=lambda x: os.path.getsize(x), reverse=True)
			
	#Create a set of reference blast databases from the input genome chunks
	def make_reference_blast_databases(self):		
		self.ref_blast_databases = []
		with multiprocessing.Pool(self.threads) as pool:
			for r in pool.imap_unordered(make_one_db, self.reference_genome_chunks):	
				self.ref_blast_databases.append(r)
		
		#Sort by file size, descending; useful for load balancing as best as possible
		self.ref_blast_databases.sort(key=lambda x: os.path.getsize(f'{x}.nsq'), reverse=True)
		
	#This code will have to be aware of genome split boundaries so it can encode which sub-file sequences ought to be retrieved from
	#For short sequences, this means just knowing which file each short sequence is contained within / need a record of: 
	#seqid : is_long + short_seq : my_file + long_seq : break starts ends
	def ref_blast(self, queries, targets):
		args = []
		for query_file in queries:
			for blast_db in targets:
				next_arg = (query_file, blast_db,)
				args.append(next_arg)
		
		mod1_gff = os.path.join(self.wd, 'module1', 'Module1_homology_hits_against_genome_gff.txt')
		mod1_gff_filt= os.path.join(self.wd, 'module1', 'Module1_homology_hits_against_genome_gff_filtered.txt')
		mod1_fasta = os.path.join(self.wd, 'module1', 'Module1_homology_hits_against_genome_fa.txt')
		mod1_fasta_filt = os.path.join(self.wd, 'module1', 'Module1_homology_hits_against_genome_fa_filtered.txt')
		json_outfile = os.path.join(self.wd, 'module1', 'Module1_homology_hits_against_genome_json.txt')
		
		o1 = open(mod1_fasta, 'w')
		o2 = open(mod1_gff, 'w')
		o3 = open(mod1_fasta_filt, 'w')
		o4 = open(mod1_gff_filt, 'w')
		
		#could use a progress tracker.
		
		#We put a bunch if info into this function and then filter what we need
		final_json = {}
		with multiprocessing.Pool(self.threads) as pool:
			for og_seqid, this_result, fasta, gff3, keeps in pool.imap_unordered(blast_ref_vs_genome, args):
				for seqid in keeps:
					this_result[seqid]['sequence_kept_after_overlaps'] = keeps[seqid].astype(int).tolist()
					for i in range(0, keeps[seqid].shape[0]):
						print(fasta[seqid][i], file = o1)
						print(gff3[seqid][i], file = o2)
						if keeps[seqid][i]:
							print(fasta[seqid][i], file = o3)
							print(gff3[seqid][i], file = o4)

					
				final_json[og_seqid] = this_result
				
		o1.close()
		o2.close()
		o3.close()
		o4.close()
		
		with open(json_outfile, 'w', encoding = 'ascii') as out:
			json.dump(final_json, out, indent = 4)
			
		return json_outfile, mod1_fasta, mod1_fasta_filt, mod1_gff, mod1_gff_filt
	
	
	#Functionally module 1
	def genome_homology(self):
		checkf = os.path.join(self.wd, 'checkpoints', 'Module1_homology_hits_against_genome_json.txt')
		checkf2 = os.path.join(self.wd, 'checkpoints', 'Module1_homology_hits_against_genome_fa.txt')
		checkf3 = os.path.join(self.wd, 'checkpoints', 'Module1_homology_hits_against_genome_fa_filtered.txt')
		checkf4 = os.path.join(self.wd, 'checkpoints', 'Module1_homology_hits_against_genome_gff.txt')
		checkf5 = os.path.join(self.wd, 'checkpoints', 'Module1_homology_hits_against_genome_gff_filtered.txt')
		
		if not os.path.exists(checkf):
			print("BLASTing reference sequences against genome...")
			
			#Split the genome into #threads chunks
			self.reference_genome_chunks = self.gs.approx_even()
			self.make_reference_blast_databases()
			
			json_outfile, mod1_fasta, mod1_fasta_filt, mod1_gff, mod1_gff_filt = self.ref_blast(self.species_libraries, 
																								self.ref_blast_databases)
			
			#Compress beforehand to make copy lighter weight
			compress(json_outfile, self.threads)
			if os.path.exists(f'{json_outfile}.gz'):
				shutil.copy(f'{json_outfile}.gz', f'{checkf}.gz')
			else:
				shutil.copy(json_outfile, checkf)
			
			shutil.copy(mod1_fasta, checkf2)
			shutil.copy(mod1_fasta_filt, checkf3)
			shutil.copy(mod1_gff, checkf4)
			shutil.copy(mod1_gff_filt, checkf5)
			
		else:
			print('Reflib vs. genome homology search already completed.')
			
		self.ref_homolog_file = checkf
		
	#Prep work for module 2, assumes GRF has already been run
	def make_library_databases(self):
		self.library_blast_dbs = []
	
		self.get_species_ref()
				
		ok_threads = min([self.threads, len(self.species_libraries)])
		with multiprocessing.Pool(ok_threads) as pool:
			for r in pool.imap_unordered(make_one_db, self.species_libraries):
				self.library_blast_dbs.append(r)
				
		#Sort by size descending
		self.library_blast_dbs.sort(key=lambda x: os.path.getsize(f'{x}.nsq'), reverse=True)
		
	def blast_from_json(self, json_file):
		print(f'BLASTing {json_file}')
		json_manager = json_loader(working_dir = self.wd)
		
		json_manager.load_json(json_file, get_names = True)
		
		args = [(w, self.library_blast_dbs[0],) for w in json_manager.workloads]
		json_manager.workloads = None
		
		clean_json_data = {}
		
		with multiprocessing.Pool(self.threads) as pool:
			for source_file, non_homolog_indices in pool.imap_unordered(json_blast, args):
				clean_json_data[source_file] = {}
				#original = json_manager.json_data[source_file]
				for seqid in non_homolog_indices:
					clean_json_data[source_file][seqid] = {'seq_length':json_manager.json_data[source_file][seqid]['seq_length'],
															'chunking_offset':json_manager.json_data[source_file][seqid]['chunking_offset'],
															'seq_start_incl_tsd':[],
															'seq_stop_incl_tsd':[],
															'tsd1_size':[],
															'tsd2_size':[],
															'tir1_size':[],
															'tir2_size':[]}
															
					for index in non_homolog_indices[seqid]:
						clean_json_data[source_file][seqid]['seq_start_incl_tsd'].append(json_manager.json_data[source_file][seqid]['seq_start_incl_tsd'][index])
						clean_json_data[source_file][seqid]['seq_stop_incl_tsd'].append(json_manager.json_data[source_file][seqid]['seq_stop_incl_tsd'][index])
						clean_json_data[source_file][seqid]['tsd1_size'].append(json_manager.json_data[source_file][seqid]['tsd1_size'][index])
						clean_json_data[source_file][seqid]['tsd2_size'].append(json_manager.json_data[source_file][seqid]['tsd2_size'][index])
						clean_json_data[source_file][seqid]['tir1_size'].append(json_manager.json_data[source_file][seqid]['tir1_size'][index])
						clean_json_data[source_file][seqid]['tir2_size'].append(json_manager.json_data[source_file][seqid]['tir2_size'][index])
					
					#free up some space
					json_manager.json_data[source_file][seqid] = None
					
		return clean_json_data	
		
	def blast_de_novo(self):
		de_novo_directory = os.path.join(self.wd, 'checkpoints')
		grf_json     = os.path.join(de_novo_directory, 'GRF_json.txt')
		tirvish_json = os.path.join(de_novo_directory, 'TIRVish_json.txt')
		
		has_grf = os.path.exists(grf_json)
		has_tirv = os.path.exists(tirvish_json)
		
		if has_tirv:
			checkf = os.path.join(self.wd, 'checkpoints', 'TIRVish_json_no_homologs.txt')
			
			if not os.path.exists(checkf):
				clean_json = self.blast_from_json(tirvish_json)
				non_homolog_tirvish_json_file = os.path.join(self.wd, 'module1', 'TIRVish_json_no_homologs.txt')
				with open(non_homolog_tirvish_json_file, 'w', encoding='ascii') as out:
					json.dump(clean_json, out, indent = 4)
					
				shutil.copy(non_homolog_tirvish_json_file, checkf)
				compress(checkf, self.threads)
			else:
				print('TIRVish homology search already completed.')
			
			self.tirvish_homologs = checkf
			
		if has_grf:
			checkf = os.path.join(self.wd, 'checkpoints', 'GRF_json_no_homologs.txt')
			if not os.path.exists(checkf):
				clean_json = self.blast_from_json(grf_json)
				non_homolog_grf_json_file = os.path.join(self.wd, 'module1', 'GRF_json_no_homologs.txt')
				with open(non_homolog_grf_json_file, 'w', encoding='ascii') as out:
					json.dump(clean_json, out, indent = 4)
					
				shutil.copy(non_homolog_grf_json_file, checkf)
				compress(checkf, self.threads)
			else:
				print('GRF homology search already completed.')
			
			self.grf_homologs = checkf
			
	def de_novo_homology(self):
		self.make_library_databases()
		self.blast_de_novo()
		
