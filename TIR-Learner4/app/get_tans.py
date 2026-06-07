import pyfastx
import numpy as np

#This needs to handle multiple sequences in the case of the short chunks.
#ista and isn just need to be dicts and seqid need to be a part of what's passed
class tan_worker:
	def __init__(self, genome_chunk, keep_sequences = False):
		self.genome_file = genome_chunk
		self.is_ta = None
		self.is_n = None
		
		self.keep = keep_sequences
		self.seq_dict = None
		self.seqlens = None

		self.prep()

	def prep(self):
		if self.keep:
			self.seq_dict = {}
			
		self.is_ta = {}
		self.is_n = {}
		self.seqlens = {}
			
		for record_tuple in pyfastx.Fasta(self.genome_file, build_index = False):
			name, seq = record_tuple[0], record_tuple[1]
			self.seqlens[name] = len(seq)
			enc = np.frombuffer(seq.encode(encoding = 'ascii'), dtype=np.int8)
			if self.keep:
				self.seq_dict[name] = seq
			seq = None
			
			self.is_ta[name] = np.zeros(enc.shape[0]+1, dtype = np.int32)
			self.is_n[name] = np.zeros(enc.shape[0]+1, dtype = np.int32)
			self.is_ta[name][1:] = np.cumsum(np.logical_or(enc == 65, enc == 84), dtype = np.int32)
			self.is_n[name][1:] = np.cumsum(enc == 78, dtype = np.int32)
		
	def check_acceptable_tans(self, seqid, starts, stops, min_seqlen = 0, max_ta_pct = 1, max_n_pct = 1):
		seqlens = stops - starts + 1
		acceptable_lengths = seqlens > min_seqlen
		acceptable_ta = (self.is_ta[seqid][stops] - self.is_ta[seqid][starts-1]) < seqlens * max_ta_pct
		acceptable_n  = (self.is_n[seqid][stops]  - self.is_n[seqid][starts-1])  < seqlens * max_n_pct
		
		passing = np.logical_and.reduce([acceptable_lengths, acceptable_ta, acceptable_n])
		
		return passing