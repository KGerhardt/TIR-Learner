#!/usr/app/env python3
# -*- coding: utf-8 -*-
# Kenji Gerhardt (kenji.gerhardt@gmail.com)
# Tianyu (Sky) Lu (tianyu@lu.fm)
# Dec. 8, 2025

import argparse
import os

VERSION = "v1.0"
INFO = "by Kenji Gerhardt, released under GPLv3"

#Options
def options():
	parser = argparse.ArgumentParser(prog="TIR-Learner",
									 description="TIR-Learner is an ensemble pipeline for Terminal Inverted Repeat "
												 "(TIR) transposable elements annotation in eukaryotic genomes")
												 
	parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION} {INFO}")

	parser.add_argument("-g", "-f", "--genome_file", help="Genome file in fasta format",
						type=str, required=True)

	parser.add_argument("-s", "--species", default = None, help="Check for homologous TIRs based on one of the following: \
																\"maize\", \"rice\". Pass \"others\" to disable homology check.")

	parser.add_argument("-l", "--length", help="Max length of TIR (Optional)", type=int, default=5000)

	parser.add_argument("-p", "-t", "--cpu", "--processors",
						help="Number of parallel processes. Default 1.", type=int, default=1, dest='processors')

	parser.add_argument("-k", "--cnn_ratio",
						help="CPU cores (OMP threads) per CNN process. Number of concurrent CNN "
							 "processes = -p // k (min 1), and each holds its own model copy, so larger "
							 "k means fewer model copies / lower peak RAM. Min 2. Default 4.",
						type=int, default=4, dest='cnn_ratio')

	parser.add_argument("-o", "--directory", "--output_dir", help=f"Where to place results. \
						Also used as the working directory. \
						Will be created if it does not exist. \
						Default location is {os.path.join(os.getcwd(), 'TIR_Learner_working_directory')}.",
						type=str, default='TIR_Learner_working_directory')

	parser.add_argument('--skip_tirvish', action = 'store_true',
						help = 'Skip running TIRvish or omit existing TIRvish results from post-processing.')
	parser.add_argument('--skip_grf', action = 'store_true',
						help = 'Skip running GRF or omit existing GRF results from post-processing.')

	parser.add_argument('--existing_tirvish', default = None,
						help = 'Supply an existing TIRvish output JSON from a previous TIR-Learner run.')
	parser.add_argument('--existing_grf', default = None,
						help = 'Supply an existing GRF output JSON from a previous TIR-Learner run.')

	# Deprecated v3 CLI flags retained for backwards compatibility (e.g. EDTA wrapper).
	# Accepted but unused; -w/--working_dir is treated as an alias of -o/--directory if -o is not given.
	parser.add_argument("-n", "--genome_name", type=str, default=None, help=argparse.SUPPRESS)
	parser.add_argument("-m", "--mode", type=str, default=None, help=argparse.SUPPRESS)
	parser.add_argument("-w", "--working_dir", type=str, default=None, help=argparse.SUPPRESS)
	parser.add_argument("-c", "--checkpoint_dir", nargs='?', const="auto", default=None, help=argparse.SUPPRESS)
	parser.add_argument("--verbose", action="store_true", help=argparse.SUPPRESS)
	parser.add_argument("-d", "--debug", action="store_true", help=argparse.SUPPRESS)
	parser.add_argument("--grf_path", type=str, default=None, help=argparse.SUPPRESS)
	parser.add_argument("--gt_path", type=str, default=None, help=argparse.SUPPRESS)
	parser.add_argument("-a", "--additional_args", type=str, nargs="+", default=None, help=argparse.SUPPRESS)

	# see prog_const for what additional args are acceptable

	parsed_args = parser.parse_args()

	return parsed_args

#Manager to handle option parsing and invoke the rest of the program. The true main behavior of TIR-Learner is in app/main.py
def main():	
	parsed_args = options()
	
	ok_species = ['rice', 'maize']
	genome_file = parsed_args.genome_file

	# -w/--working_dir is a v3 alias; only honor it if -o/--directory was left at default
	directory = parsed_args.directory
	if parsed_args.working_dir is not None and directory == 'TIR_Learner_working_directory':
		directory = parsed_args.working_dir

	threads = parsed_args.processors
	species = parsed_args.species

	# Normalize species: EDTA passes "Rice"/"Maize"/"others"; v3 also accepted "others" to skip homology
	if species is not None:
		species = species.lower()
		if species == 'others':
			species = None

	tir_max_length = parsed_args.length

	ok_to_continue = True

	if species is not None:
		if species not in ok_species:
			print(f'Your supplied species {species} is not in the acceptable species list, which is:')
			for s in ok_species:
				print(f'\t{s}')
			ok_to_continue = False
						
	# CNN cores-per-process. Drives both OMP_NUM_THREADS and the CNN pool size (threads // cnn_ratio).
	# Clamped to >= 2 so it can never collapse to a single OMP thread on a multi-process run.
	cnn_ratio = max(2, parsed_args.cnn_ratio)

	if ok_to_continue:
		if threads == 1:
			os.environ['OMP_NUM_THREADS'] = '1'
		else:
			#We handle most of the parallelization in python; give each CNN process cnn_ratio OMP threads
			os.environ['OMP_NUM_THREADS'] = str(cnn_ratio)
		os.environ["KERAS_BACKEND"] = "torch"  # use pytorch as keras backend
		os.environ["KMP_WARNINGS"] = '0'  # mute all OpenMP warnings
		
		os.environ["CUDA_VISIBLE_DEVICES"] = "," #force torch to use CPU, not GPU

	
		print('Loading program resources. This usually takes about a minute.')
		from app.main import newTL
		print('Program resources loaded.')
		print('')
	
		skip_t = parsed_args.skip_tirvish
		skip_g = parsed_args.skip_grf
		exist_t = parsed_args.existing_tirvish
		exist_g = parsed_args.existing_grf
	
		TIRLearner_instance = newTL(
								genome_file_path = genome_file,
								TIR_length = tir_max_length,
								processors = threads,
								cnn_ratio = cnn_ratio,
								species = species,
								wd = directory,
								skip_tirvish = skip_t,
								skip_grf = skip_g, 
								existing_tirvish = exist_t,
								existing_grf = exist_g
								)
		TIRLearner_instance.run()

if __name__ == "__main__":
	main()
