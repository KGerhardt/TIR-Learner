#!/usr/app/env python3
# -*- coding: utf-8 -*-
# Kenji Gerhardt
# Tianyu (Sky) Lu (tianyu@lu.fm)
# 2025-12-02

import sys
import os

from .genomeSplitter import genomeSplitter

from .tirvish_new import TIRvish_manager
from .grf_new import GRF_manager
from .cnn_new import CNN_manager
from .blast_new import blaster

from .output_compressor import compress, decompress

import shutil


'''
Notes for future development:

The overwhelmingly likely use-case for this code is GRF + TIRVish + CNN with no homology search. There is
still some small but meaningful speedup to be gained:

Each genome chunk is already loaded twice for TIRvish and GRF as part of the tan_checker code
The GRF sequences are loaded once as part of the GRF output filtering stage and once as part of the CNN processing step
TIRVish has some rare cases where its performance is terrible on just a few genome chunks and takes 10s of minutes to process 5Mbp

Running GRF+TIRvish concurrently and immediately applying CNN to GRF results to avoid double I/O is probably best-case. However, this is likely
a bunch of work for *very* small time gains in most cases.
'''

class newTL:
	def __init__(self, genome_file_path: str, TIR_length: int = 5_000,
				 processors: int = 1, species = None, wd = 'TIR_Learner_working_directory',
				 extension_size = 20, chunk_size = 5_000_000, olap = 7_500, skip_tirvish = False, 
				 skip_grf = False, existing_tirvish = None, existing_grf = None):
		
		self.threads = processors
		
		self.gf = genome_file_path
		self.gfi = f'{self.gf}.fxi'
		self.genome_name = os.path.basename(self.gf)
		
		self.wd = wd
		self.split_dir = os.path.join(self.wd, 'split_genome')
		self.blast_dir = os.path.join(self.wd, 'module1')
		self.check_dir = os.path.join(self.wd, 'checkpoints')
		self.active_dir = os.path.join(self.wd, 'current_results')
		self.final_dir = os.path.join(self.wd, 'TIR-Learner-Result')
		
		self.chunk = chunk_size
		
		self.extend = extension_size
		
		self.TIR_length = TIR_length
		self.olap = self.TIR_length + extension_size
		
		self.species = species

		self.gsplit = genomeSplitter(genome_file = self.gf,
									output_directory = self.split_dir,
									chunk_size = chunk_size,
									#Use a safe overlap size, Ensure that extension is possible later
									overlap_size = self.olap,
									minimum_seq_size = self.TIR_length + 500,
									procs = self.threads,
									smart = False,
									post_index = False,
									verbose = True,
									overwrite = False,
									do_bedtools_prep = True)
		
		self.existing_grf = existing_grf
		self.existing_tirvish = existing_tirvish
		
		self.skip_tirvish = skip_tirvish
		self.skip_grf = skip_grf
									
		self.split_genome_files = None
		
		self.seqlens = None
		
		self.tirvish_file = None
		self.grf_file = None
		
		self.blast_file = None
		
		self.cnn_file = None
		
	def prepare_directory(self, this_dir):
		if not os.path.exists(this_dir):
			os.makedirs(this_dir, exist_ok = True)
			
	def dir_prep_pre(self):
		self.prepare_directory(self.split_dir)
		self.prepare_directory(self.check_dir)
		self.prepare_directory(self.active_dir)
		
	#Execute genome splitter if needed
	def scan_and_split_genome(self):
		self.split_genome_files = self.gsplit.run()
		self.seqlens = self.gsplit.seqs_and_lengths
		
	#Runs GRF, processes results into a JSON file
	def GRF(self):
		if self.existing_grf is None:
			self.grf_file = GRF_manager(self.split_genome_files, self.seqlens, self.active_dir, self.check_dir, self.olap, self.chunk, threads = self.threads)
			compress(self.grf_file, threads = self.threads)
		else:
			if os.path.exists(self.existing_grf):
				recheck = os.path.join(self.check_dir, os.path.basename(self.existing_grf))
				shutil.copy(self.existing_grf, recheck)
				self.grf_file = decompress(recheck)
			else:
				print(f'Supplied GRF JSON {self.existing_grf} could not be found.')
				
					
	#Runs TIRvish, processes results into a json file
	def TIRvish(self):
		if self.existing_tirvish is None:
			self.tirvish_file = TIRvish_manager(self.split_genome_files, self.seqlens, self.active_dir, self.check_dir, self.olap, self.chunk, threads = self.threads)
			compress(self.tirvish_file, threads = self.threads)
		else:
			if os.path.exists(self.existing_tirvish):
				recheck = os.path.join(self.check_dir, os.path.basename(self.existing_tirvish))
				shutil.copy(self.existing_tirvish, recheck)
				self.tirvish_file = decompress(recheck)
			else:
				print(f'Supplied tirvish JSON {self.existing_tirvish} could not be found.')
	
	#This is 95% implemented
	def blast(self):
		#Nice, simple way to make this check; we check for OK species in main
		if self.species is not None:
			blast_manager = blaster(reference_genome = self.gf, 
									species = self.species, 
									working_dir = self.wd, 
									threads = self.threads)
									
			blast_manager.genome_homology()
			blast_manager.de_novo_homology()

	
	def CNN(self):
		#Use half threads and give each CNN process 2 OMP threads; balances speed and memory
		manager = CNN_manager(self.tirvish_file, 
					self.grf_file, 
					self.wd, 
					max([self.threads // 2, 1]))
		
		final_fa, final_g3, final_fa_filt, final_g3_filt, final_tirvish, final_grf = manager.run()
		
		self.prepare_directory(self.final_dir)
		#Move all TIRVish materials
		if self.tirvish_file is not None:
			for f in os.listdir(self.check_dir):
				if 'TIRVish_json' in f and not f.endswith('.gz'):
					compress(os.path.join(self.check_dir, f), threads = self.threads)
				
			base_file = os.path.join(self.check_dir, 'TIRVish_json.txt.gz')
			homolog_file = os.path.join(self.check_dir, 'TIRVish_json_no_homologs.txt.gz')
			base_dest = os.path.join(self.final_dir, 'TIRVish_json.txt.gz')
			homolog_dest = os.path.join(self.final_dir, 'TIRVish_json_no_homologs.txt.gz')
			
			shutil.move(base_file, base_dest)
			if os.path.exists(homolog_file):
				shutil.move(homolog_file, homolog_dest)
				
			#compress(final_tirvish, threads = self.threads)
			cnn_dest = os.path.join(self.final_dir, 'post_CNN_TIRVish_json.txt.gz')
			shutil.move(f'{final_tirvish}.gz', cnn_dest)
		
		#Move all GRF materials	
		if self.grf_file is not None:
			for f in os.listdir(self.check_dir):
				if 'GRF_json' in f and not f.endswith('.gz'):
					compress(os.path.join(self.check_dir, f), threads = self.threads)
					
			base_file = os.path.join(self.check_dir, 'GRF_json.txt.gz')
			homolog_file = os.path.join(self.check_dir, 'GRF_json_no_homologs.txt.gz')
			base_dest = os.path.join(self.final_dir, 'GRF_json.txt.gz')
			homolog_dest = os.path.join(self.final_dir, 'GRF_json_no_homologs.txt.gz')
			
			shutil.move(base_file, base_dest)
			if os.path.exists(homolog_file):					
				shutil.move(homolog_file, homolog_dest)
			
			#compress(final_grf, threads = self.threads)
			cnn_dest = os.path.join(self.final_dir, 'post_CNN_GRF_json.txt.gz')
			shutil.move(f'{final_grf}.gz', cnn_dest)

		#Move homology results if needed
		for f in os.listdir(self.check_dir):
			if f.startswith('Module1_homology_hits_against_genome_'):
				shutil.move(os.path.join(self.check_dir, f), os.path.join(self.final_dir, f))
		
		#Move all finalized CNN results
		shutil.move(final_fa, os.path.join(self.final_dir, os.path.basename(final_fa)))
		shutil.move(final_fa_filt, os.path.join(self.final_dir, os.path.basename(final_fa_filt)))
		shutil.move(final_g3, os.path.join(self.final_dir, os.path.basename(final_g3)))
		shutil.move(final_g3_filt, os.path.join(self.final_dir, os.path.basename(final_g3_filt)))
		
	def clean_up(self):
		print('')
		print('Cleaning up temporary directories')
		
		print('Removing partial results directory')
		shutil.rmtree(self.active_dir)
		print('Removing genome chunk directory')
		shutil.rmtree(self.split_dir)
		if os.path.exists(self.blast_dir):
			print('Removing genome homology directory')
			shutil.rmtree(self.blast_dir)
			
		print(f'Removing checkpoints directory (all important checkpoints are still kept in {self.final_dir})')
		shutil.rmtree(self.check_dir)
		
	def run(self):
		#Prep work - always run
		self.dir_prep_pre()

		#This is a prereq for tirvish/grf, so will always be run
		self.scan_and_split_genome()

		#These will be run irrespective of module 1-3 or 4
		if not self.skip_tirvish:
			self.TIRvish()
		if not self.skip_grf:
			self.GRF()
		if self.species is not None:
			self.blast()

		self.CNN()

		self.clean_up()

