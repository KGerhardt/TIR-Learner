import os
import multiprocessing
import pyfastx
import sqlite3
import argparse

#import shutil
import subprocess

#import numpy as np

def options():
	parser = argparse.ArgumentParser(
		description='''Tool and Python API for splitting a genome into chunks.
		
	* Long sequences (larger than --chunk_size bp) are split into the fewest chunks so 
	  no chunk is larger than --chunk_size bp including an --overlap_size bp overlap.
	* Short sequences are grouped into multiFASTA files around --chunk_size bp - may go very slightly over.''',
		formatter_class=argparse.RawDescriptionHelpFormatter
	)
	
	# Required arguments
	parser.add_argument('--genome', required=True,
						help='Path to genome FASTA file.')
	parser.add_argument('--output_directory', required=True,
						help='Output directory path for split file results.')
	
	# Integer arguments with validation
	parser.add_argument('--chunk_size', default = 250_000_000, type=int,
						help='Size of processing chunks in base pairs. Default 250 million bp.')
	parser.add_argument('--overlap_size', default = 1_000_000, type=int,
						help='Size of overlap between chunks in base pairs. Default 1 million bp')
	parser.add_argument('--processors', type=int, default=1,
						help='Number of processors to use')
	
	# Boolean flags
	parser.add_argument('--verbose', action='store_true',
						help='Enable verbose output')
	parser.add_argument('--index_outputs', action='store_true',
						help='Create Pyfastx indices of split genome files to facilitate later processing.')
	parser.add_argument('--smart', action='store_true',
						help='Make an effort to divide the genome into --processors chunks so that later parallel processing can be ~1 process/chunk. No guarantees.')
	parser.add_argument('--overwrite', action='store_true',
						help='Delete existing genomeSplitter outputs in the output directory, if any are found.')


	args = parser.parse_args()

	# Add validation checks
	if args.chunk_size <= 0:
		parser.error('--chunk_size must be a positive integer')
	if args.overlap_size < 0:
		parser.error('--overlap_size must be a non-negative integer')
	if args.overlap_size >= args.chunk_size:
		parser.error('--overlap_size must be smaller than --chunk_size')
	if args.processors < 1:
		parser.error('--processors must be at least 1')
		
	return args

	
class manual_fai:
	def __init__(self):
		self.current_offset = 0
		self.records = []
	
	def add_record(self, name, unformatted_seq):
		line_size = 70
		line_bytes = line_size + 1 #bp per line + '\n'
		header_offset = len(name) + 2 #'>' + length of seqid + '\n'
		seqlen = len(unformatted_seq)
		formatted_seq = format_num_bp(unformatted_seq, line_size)
		next_record = (name, seqlen, self.current_offset + header_offset, line_size, line_bytes,)
		self.records.append(next_record)
		self.current_offset = self.current_offset + header_offset + len(formatted_seq) + 1
		
		return formatted_seq
		
	def purge(self):
		self.current_offset = 0
		self.records = []
	
class genomeSplitter:
	def __init__(self, genome_file, output_directory, chunk_size = 250_000_000, 
				overlap_size = 1_000_000, minimum_seq_size = None, procs = 1, smart = False, post_index = False, 
				verbose = False, overwrite = False, quiet = False, do_bedtools_prep = False):
		self.path = os.path.abspath(genome_file)
		self.index = f'{self.path}.fxi'
		
		self.outdir = output_directory
		
		self.seqs_and_lengths = None
		
		self.minsize = minimum_seq_size
		self.chunk = chunk_size
		self.olap = overlap_size
		
		self.short_plan = None
		self.long_plan = None
		self.overall_split_plan = None
		
		self.smart = smart
		self.threads = procs
		self.verbose = verbose
		self.quiet = quiet
		self.do_output_index = post_index
		self.overwrite = overwrite
		
		self.output_files = None
		
		self.do_bedtools_prep = do_bedtools_prep
		#self.bed_prep_dir = bedtools_prep_dir
		
	def get_seqlens(self):
		if self.minsize is not None:
			retrieve = f'SELECT chrom, slen FROM seq WHERE slen >= {self.minsize}'
		else:
			retrieve = f'SELECT chrom, slen FROM seq'
	
		conn = sqlite3.connect(self.index)
		curs = conn.cursor()
		summary = curs.execute(retrieve).fetchall()
		curs.close()
		conn.close()
		self.seqs_and_lengths = dict(summary)
				
		#for row in summary:
		#	seq_name = row[0]
		#	seq_length = row[1]
		#	self.seqs_and_lengths[seq_name] = seq_length
		
	#Create pyfastx index for the genome if it doesnt exist; load sequence lengths for planning.
	def index_and_summarize(self):
		if os.path.exists(self.index):
			if not self.quiet:
				print(f'Pyfastx index already found for {os.path.basename(self.path)}.')
				print(f'Loading sequence summary.')
		else:
			if not self.quiet:
				print(f'Beginning Pyfastx index process on {os.path.basename(self.path)}.')
				print(f'This should only take a few seconds per billion bp in your genome.')
			fa = pyfastx.Fasta(self.path, build_index = True)
			if not self.quiet:
				print('Indexing complete.')
		
		self.get_seqlens()
		
		if not self.quiet:
			print(f'Sequence summarized.')

	#Divide long sequences into even chunks with overlap as close to chunk size as possible without going over
	def evenly_chunk_long_sequences(self, longs):
		long_split_plan = {}
		if self.verbose and len(longs) > 0:
			print('')
			print(f'Long sequence splits have been planned. Here is the summary:')
			print('')
		for seqid in longs:
			split_count = 0
			long_split_plan[seqid] = []
			length = longs[seqid]
			
			num_chunks = length // self.chunk #Guesstimate of the initial number of chunks to make
			adjusted_length = length + ((num_chunks-1) * self.olap) #We effectively extend the sequence by self.olap x num chunks
			while int(adjusted_length / num_chunks) > self.chunk:
				num_chunks += 1
				adjusted_length = length + ((num_chunks-1) * self.olap)
				
			adjusted_chunk_size = int(length / num_chunks) + self.olap
			
			start, end = 0, 0 #initialize
			while end < length: #simplify loop check
				end = min(start + adjusted_chunk_size, length)
				long_split_plan[seqid].append((start, end,))
				start = end - self.olap
				split_count += 1
				
			if self.verbose:
				print(f'\tLong sequence {seqid} ({length} bp) was split into {split_count} chunks of {adjusted_chunk_size} bp')
		
		if self.verbose:
			print('')
		
		return long_split_plan
			
	#Aggregate small sequences into the fewest groups not larger than chunk size; pack each group with about same number of bp
	def aggregate_small_sequences(self, shorts):
		total_seqlen = 0
		for seqid in shorts:
			total_seqlen += shorts[seqid]
		
		#Ceiling function
		minimum_number_of_files = int(round((total_seqlen / self.chunk) + 0.5, 0))
		
		#Descending seqlen size sorted dict
		shorts = dict(sorted(shorts.items(), key=lambda item: item[1], reverse = True))
		
		short_split_plan = {}
		chunk_length_record = {}
		for i in range(0, minimum_number_of_files):
			short_split_plan[i] = []
			chunk_length_record[i] = 0
		
		shortest_chunk = 0
		for seqid in shorts:
			seqlen = shorts[seqid]
			short_split_plan[shortest_chunk].append(seqid)
			chunk_length_record[shortest_chunk] += seqlen
			shortest_chunk = min(chunk_length_record, key=chunk_length_record.get)			
		
		keyset = list(short_split_plan.keys())
		for k in keyset:
			if len(short_split_plan[k]) == 0:
				delete = short_split_plan.pop(k)
		
		if self.verbose and len(shorts) > 0:
			print('Short sequences have been aggregated. Here is the summary:')
			print('')
			for i in short_split_plan:
				print(f'\tSplit {i+1} contains {chunk_length_record[i]} bp across {len(short_split_plan[i])} sequences.')
			print('')
				
		return short_split_plan
	
	def one_per_proc(self):
		if self.verbose:
			print(f'Trying to split the genome into {self.threads} chunks, one per processor')
		total_genome_size = sum(self.seqs_and_lengths.values())
		self.chunk = int(round(((total_genome_size / self.threads) + 0.5), 0)) + self.threads
		self.olap = self.chunk // 25
		if self.verbose:
			print(f'New chunk size {self.chunk} overlap is 4% of chunk size at {self.olap} bp')
	
	#Divide a genome into exactly self.threads chunks of unidivided sequences as close to evenly as possible
	def approx_even(self):
		self.index_and_summarize()
		
		grouping_size = sum(self.seqs_and_lengths.values()) / self.threads
		self.chunk = grouping_size
		
		self.short_plan = self.aggregate_small_sequences(self.seqs_and_lengths)
		
		self.overall_split_plan = []
		
		for group in self.short_plan:
			next_group = ('short', group, self.short_plan[group], self.path, self.outdir, self.do_output_index,)
			self.overall_split_plan.append(next_group)
		
		self.long_plan = None
		self.short_plan = None
		
		if not self.quiet:
			print(f'Genome split plan complete.')
			
		self.compare_or_skip()
		
		return self.output_files
		
		
	#Assess the genome and create a plan to split it elegantly
	def prepare_split_plan(self, chunk_size = None, overlap_size = None):
		#Check manual sizes are OK.
		if chunk_size is not None or overlap_size is not None:
			if isinstance(chunk_size, int) and isinstance(overlap_size, int):
				if not self.quiet:
					print(f'Using newly supplied chunk size {chunk_size} and overlap size {overlap_size}')
				self.chunk = chunk_size
				self.olap = overlap_size
			else:
				if not self.quiet:
					print(f'It looks like you tried to supply a new size of chunk and overlap, but they were not properly formatted.')
					print(f'These must both be integers. Supplied values:')
					print(f'\tChunk size: {chunk_size}')
					print(f'\tOverlap size: {overlap_size}')
					print('')
					print(f'Resorting to program defaults of chunk size {self.chunk} and overlap size {self.olap}')
		else: #do not try to be smart if sizes were manually supplied.
			if self.smart:
				self.one_per_proc()

		long_sequences = {}
		short_sequences = {}
		
		for seqid in self.seqs_and_lengths:
			seqlen = self.seqs_and_lengths[seqid]
				
			if seqlen > self.chunk:
				long_sequences[seqid] = seqlen
			else:
				short_sequences[seqid] = seqlen
		
		self.long_plan = self.evenly_chunk_long_sequences(long_sequences)
		self.short_plan = self.aggregate_small_sequences(short_sequences)
		
		self.overall_split_plan = []
		
		#Because of how pyfastx works, it will be quicker and smarter to load each seqid in one writer process and let it go to work
		for seqid in self.long_plan:
			next_group = ('long', seqid, self.long_plan[seqid], self.path, self.outdir, self.do_output_index,)
			self.overall_split_plan.append(next_group)
		
		for group in self.short_plan:
			next_group = ('short', group, self.short_plan[group], self.path, self.outdir, self.do_output_index,)
			self.overall_split_plan.append(next_group)
		
		self.long_plan = None
		self.short_plan = None
		
		if not self.quiet:
			print(f'Genome split plan complete.')
		
	def prep_outdir(self):
		if not os.path.exists(self.outdir):
			print(f'Making directory {self.outdir}')
			os.makedirs(self.outdir, exist_ok = True)

	def execute_split(self):
		if not self.quiet:
			print('Executing genome split.')
		total_outputs = len(self.overall_split_plan)
		if self.verbose:
			print(f'Genome {os.path.basename(self.path)} will be split into {total_outputs} chunked outputs.')
			print(f'Results will be found at {os.path.abspath(self.outdir)}')
			print('')
			
		self.prep_outdir()
		
		if self.verbose:
			if total_outputs < self.threads:
				print(f'The split plan has only {total_outputs} chunks. Only this many parallel processes will be used.')
				
		#Never use more threads than there are splits.
		ok_threads = min([total_outputs, self.threads])
		
		self.output_files = []
		
		work_done = 1
		with multiprocessing.Pool(ok_threads) as pool:
			for result in pool.imap_unordered(chunk_write, self.overall_split_plan):
				self.output_files.extend(result)
				if self.verbose:
					print(f'{work_done} of {total_outputs} complete.')
				work_done += 1
				
		self.output_files.sort()
		
		'''
		#This part is now manually done
		if self.do_bedtools_prep:
			print('Creating samtools fasta indices for genome fragments...')
			with multiprocessing.Pool(ok_threads) as pool:
				for result in pool.imap_unordered(make_bedtools_idx, self.output_files):
					pass
					
			print('.fai indices created!')
		'''
		
		if not self.quiet:
			print('Genome split complete.')
				
		return self.output_files
	
	def indices_only(self):
		outputs = [os.path.join(self.outdir, os.path.basename(f)) for f in self.output_files]
		completed_indices = []
		ok_threads = min([len(outputs), self.threads])
		with multiprocessing.Pool(ok_threads) as pool:
			for result in pool.imap_unordered(index_only, outputs):
				completed_indices.append(result)
	
		completed_indices.sort()
		self.output_files = completed_indices
		
		return self.output_files
	
	def compare_log_to_plan(self, log_file):
		outputs = []
		with open(log_file) as fh:
			for line in fh:
				param = line.strip().split('\t')[1]
				outputs.append(param)
		
		plan = []
		plan.append(str(self.chunk))
		plan.append(str(self.olap))
		
		plan_and_log_are_identical = outputs[0] == plan[0] and outputs[1] == plan[1]
		#Check if this run is asking for indices
		
		just_index = self.do_output_index and outputs[2] == "False"
		
		return plan_and_log_are_identical, just_index, outputs[3:]
			
	def clean_outdir(self):
		for filename in os.listdir(self.outdir):
			file_path = os.path.join(self.outdir, filename)
			if os.path.isfile(file_path):  # Check if it's a file
				os.remove(file_path)
	
	def create_log(self, log_file):
		with open(log_file, 'w') as out:
			print(f'chunk_size\t{self.chunk}', file = out)
			print(f'overlap_size\t{self.olap}', file = out)
			print(f'indices_created\t{self.do_output_index}', file = out)
			for f in self.output_files:
				print(f'{self.do_output_index}\t{f}', file = out)

	#Check if the log file exists, indicating a successful run, and compare to the existing split plan; run the split if the log is not OK
	def compare_or_skip(self):
		already_run = False
		needs_index = False
		log_file = os.path.join(self.outdir, 'genomeSplitter.log')
		if os.path.exists(log_file):
			#Silently run this
			unchanged_plan, needs_index, previous_log = self.compare_log_to_plan(log_file)
			if not self.quiet:
				print('')
				print(f'It looks like genomeSplitter was run previously in {self.outdir}')
			if self.overwrite:
				if unchanged_plan:
					if not self.quiet:
						print('Even though --overwrite was supplied, the current plan is identical to the previous run.')
						print('')
						print('Previous and current run shared parameters:')
						print(f'\tChunk size = {self.chunk}')
						print(f'\tOverlap size = {self.olap}')
						print('')
					if needs_index:
						if not self.quiet:
							print('However, this run requested Pyfastx indices for the outputs and those were not previously created.')
							print('Those will be created now.')
					else:
						if not self.quiet:
							print('There is no point in rerunning genomeSplitter unless the parameters are different.')
					already_run = True
				else:
					if not self.quiet:
						print('Different parameters were supplied to the current run. Removing old files and splitting the genome.')
					self.clean_outdir()
			else:
				if not self.quiet:
					print('')
					print(f'--overwrite was not supplied. genomeSplitter will not try to re-split this genome.')
				already_run = True
				if needs_index:
					if not self.quiet:
						print('However, this run requested Pyfastx indices for the outputs and those were not previously created.')
						print('Those will be created now.')
			
			if not self.quiet:		
				print('')
		
		if not already_run:	
			self.execute_split()
			self.create_log(log_file)
			
		else:
			self.output_files = previous_log
			self.output_files.sort()
			if needs_index:	
				self.indices_only()
				self.create_log(log_file)
		
		
	def run(self):
		self.index_and_summarize()
		self.prepare_split_plan()
		
		self.compare_or_skip()
		
		return self.output_files
		
def format_num_bp(string, size = 70):
	politely_formatted = '\n'.join([(string[i:i+size]) for i in range(0, len(string), size)])
	return politely_formatted
			
def chunk_write(args):
	short_or_long = args[0]
	genome_file = args[3]
	output_dir = args[4]
	do_index = args[5]
	#do_bed = args[6]
	#bed_dir = args[7]
	
	idx_worker = manual_fai()
	
	fa = pyfastx.Fasta(genome_file)
	
	output_files = []
	
	if short_or_long == "short":
		chunk_id = args[1]
		sequences = args[2]
		output_file = os.path.join(output_dir, f'short_chunk_{chunk_id}_offset_0.fasta')
		#my_fai = f'{output_file}.fai'
		
		with open(output_file, 'w') as out:
			for seqid in sequences:
				sequence = fa[seqid].seq.upper()

				new_id = f'{seqid};;0'
				
				writeout = idx_worker.add_record(new_id, sequence)
				print(f'>{new_id}', file = out)
				print(f'{writeout}', file = out)

		#Create fasta index
		#with open(my_fai, 'w') as out:
		#	for record in idx_worker.records:
		#		print(*record, sep = '\t', file = out)
				
		output_files.append(output_file)
				
		if do_index:
			fa = pyfastx.Fasta(output_file, build_index = True)
	else:
		seqid = args[1]
		print(f'Wrangling long sequence {seqid}')
		long_seq = fa[seqid].seq.upper()
		
		print(f'Long sequence {seqid} total bp is {len(long_seq)} and has {len(args[2])} chunks')
		
		for sequence_interval in args[2]:
			start = sequence_interval[0]
			end = sequence_interval[1]
			
			output_file = os.path.join(output_dir, f'long_chunk_{seqid}_offset_{start}.fasta')
			#my_fai = os.path.join(f'{output_file}.fai')
			new_id = f'{seqid};;{start}'
			subseq = long_seq[start:end]
			writeout = idx_worker.add_record(new_id, subseq)
			with open(output_file, 'w') as out:
				print(f'>{new_id}', file = out)
				print(f'{writeout}', file = out)
			
			#Create fasta index
			#with open(my_fai, 'w') as out:
			#	for record in idx_worker.records:
			#		print(*record, sep = '\t', file = out)
				
			#Reset the indexer
			idx_worker.purge()
			
			if do_index:
				fa = pyfastx.Fasta(output_file, build_index = True)

			output_files.append(output_file)

	return output_files

#Obsolete, replaced with a pure python equivalent of the same
def make_bedtools_idx(file):
	comm = f'samtools faidx {file} --threads 1'
	comm = comm.split()
	subprocess.run(comm)
	
	os.system(f'diff {file}.fai {file}.kenji.fai')
	
	os.remove(f'{file}.fai')
	
def index_only(filepath):
	pyfastx.Fasta(filepath, build_index = True)
	return filepath

def main():
	#Used in parallel processes
	opts = options()
	genome_file = opts.genome
	output_dir = opts.output_directory
	index_outputs = opts.index_outputs
	
	p = opts.processors
	c = opts.chunk_size
	o = opts.overlap_size
	v = opts.verbose
	s = opts.smart
	ow = opts.overwrite
	
	mn = genomeSplitter(genome_file = genome_file, 
				output_directory = output_dir, 
				chunk_size = c,
				overlap_size = o,
				procs = p,
				smart = s,
				post_index = index_outputs,
				verbose = v,			
				overwrite = ow)
				
	
				
	output_files = mn.run()
	
if __name__ == '__main__':
	main()