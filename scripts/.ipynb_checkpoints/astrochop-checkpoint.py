import os
import sys

import subprocess as sp
import numpy as np

# DEFINE METHODS AND THEIR PATH

methods = {'chainsaw': 'path_to_chainsaw'}

# DEFINE MAIN FUNCTIONS

# 1. main functions

def chop_domains(inpdb, methods=methods):
	
	predictions = []
	for method in methods:
		domains = run_predictor(inpdb, method, path = methods[method])
		predictions.append(domains)

		print(domains)
	
	domains, consensus = get_consensus(predictions)

	print(consensus)
	
	return domains, len(consensus)

def get_consensus(predictions):
	
	print('Computing consensus')
	
	predictions = np.array(predictions).T
	consensus = [None for i in predictions]
	
	for i, residue in enumerate(predictions):
		residue = [j for j in residue if j is not None]
		if len(residue) > 0:
			score = np.median(residue)
			if score > 0:
				score -= 0.1
		else:
			score = None

		try:
			consensus[i] = int(round(score, 0))
		except:
			consensus[i] = score
	
	domains = convert_map_to_interval(consensus)				
	print(' ... Found {} domains: {}'.format(len(domains), domains))
	
	return domains, consensus

# 2. predictors

def run_predictor(inpdb, method, path):
	
	def run_merizo(inpdb, path):
		
		print('Running Merizo')
		run_merizo = sp.Popen(['python3', path, '-i', inpdb], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
		stdout, stderr = run_merizo.communicate()
		
		domains = stdout.decode().strip().split('\t')[-1].split(',')
		nres = int(stdout.decode().strip().split('\t')[1])
		
		return domains, nres
	
	def run_chainsaw(inpdb, path):
		
		print('Running ChainSaw')
		outfile = '{}.chainsaw'.format(inpdb)
		
		if not os.path.isfile(outfile):
			run_chainsaw = sp.Popen(['python3', path, '--structure_file', inpdb, '--output', outfile], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
			stdout, stderr = run_chainsaw.communicate()
		
		if os.path.isfile(outfile):
			with open(outfile, 'r') as outf:
				for line in outf:
					if not line.startswith('chain_id'):
						domains = line.strip().split('\t')[-2].split(',')
						nres = int(line.strip().split('\t')[2])
						return domains, nres
		else:
			raise Exception(stderr.decode())
	
	#----------------------
	if method == 'merizo':
		domains, nres = run_merizo(inpdb, path)
	elif method == 'chainsaw':
		domains, nres = run_chainsaw(inpdb, path)
	
	print(' ... Found {} domains: {}'.format(len(domains), domains))
	
	resmap = convert_interval_to_map(domains, nres)
	
	return resmap

# 3. helper functions
def convert_map_to_interval(resmap):
	
	domains = {}
	for i, pred in enumerate(resmap):
		if pred is not None:
			if pred not in domains:
				domains[pred] = []
			domains[pred].append(i)
	
	for domain in domains:
		parts = [[]]
		for i in domains[domain]:
			if len(parts[-1]) == 0 or i - parts[-1][-1] == 1:
				parts[-1].append(i)

			else:
				parts.append([i])

		for i in range(len(parts)):
			parts[i] = '{}-{}'.format(min(parts[i])+1, max(parts[i])+1)

		if len(parts) > 1:
			domains[domain] = '_'.join(parts)
		else:
			domains[domain] = parts[0]
	
	return list(domains.values())

def convert_interval_to_map(domains, nres):
	
	resmap = [None for i in range(nres)]
	for indx, domain in enumerate(domains):
		parts = domain.split('_')
		for part in parts:
			try:
				start, end = [int(i) for i in part.split('-')]
				if start == 0:
					start = 1
				for i in range(start-1, end):
					resmap[i] = indx+1
			except:
				pass
	
	return resmap

###########

if __name__ ==  '__main__':

	inpdb = sys.argv[1]

	domains, nres = chop_domains(inpdb)

	print('\nCOLMNS:\tinpdb\tnres\tdomains')
	print('RESULT:\t{}\t{}\t{}'.format(inpdb, nres, ','.join(domains)))