import subprocess
import os
import json

def write_formatted_json(obj, output_path):
	"""Serialize `obj` to `output_path` as compact JSON (no whitespace).

	Compact output is ~5x smaller on disk than indent=4 at identical serialization
	speed (the indentation is pure bloat for these machine-read checkpoints).
	Centralized here so the JSON format/serializer can be changed in one place.
	"""
	with open(output_path, 'w', encoding = 'ascii') as out:
		json.dump(obj, out, separators = (',', ':'))

def read_json(input_path):
	"""Load a JSON file. Counterpart to write_formatted_json that centralizes JSON
	deserialization so the parser can be swapped in one place. Compact JSON parses
	identically to indented JSON (separators only affect writing), so this reads
	files written either way.
	"""
	with open(input_path, 'r', encoding = 'ascii') as inf:
		return json.load(inf)

def compress(input_file, threads = 1):
	if os.path.exists(f'{input_file}.gz'):
		print(f'{input_file} gzip was already found. Nothing to do.')
	else:
		try:
			print(f'Attempting to compress {input_file} with pigz...')
			proc = subprocess.call(['pigz', '-6', '-k', '-p', str(threads), input_file])
		except Exception:
			print('pigz compressor not found. Defaulting to gzip.')
			try:
				print(f'Attempting to compress {input_file} with gzip...')
				proc = subprocess.call(['gzip', '-k', '-6', input_file])
			except Exception:
				print('gzip compressor not found. The file will be left uncompressed.')

def decompress(input_file, threads = 1):
	if input_file.endswith('.gz'):
		no_zip = input_file[:-3]
		if not os.path.exists(no_zip):
			try:
				subprocess.call(['pigz', '-d', '-k', '-p', str(threads), input_file])
			except Exception:
				print('pigz decompress failed')
				try:
					subprocess.call(['gunzip', '-k', input_file])
				except Exception:
					print('gunzip decompress failed.')
	else:
		print(f'File {input_file} was already unzipped')
		no_zip = input_file
		
	return no_zip