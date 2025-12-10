import subprocess
import os

def compress(input_file, threads = 1):
	if os.path.exists(f'{input_file}.gz'):
		print(f'{input_file} gzip was already found. Nothing to do.')
	else:
		try:
			print(f'Attempting to compress {input_file} with pigz...')
			proc = subprocess.call(['pigz', '-6', '-k', '-p', str(threads), input_file])
		except:
			print('pigz compressor not found. Defaulting to gzip.')
			try:
				print(f'Attempting to compress {input_file} with gzip...')
				proc = subprocess.call(['gzip', '-k', '-6', input_file])
			except:
				print('gzip compressor not found. The file will be left uncompressed.')

def decompress(input_file, threads = 1):
	if input_file.endswith('.gz'):
		no_zip = input_file[:-3]
		if not os.path.exists(no_zip):
			try:
				subprocess.call(['pigz', '-d', '-k', '-p', str(threads), input_file])
			except:
				print('pigz decompress failed')
				try:
					subprocess.call(['gunzip', '-k', input_file])
				except:
					print('gunzip decompress failed.')
	else:
		print(f'File {input_file} was already unzipped')
		no_zip = input_file
		
	return no_zip