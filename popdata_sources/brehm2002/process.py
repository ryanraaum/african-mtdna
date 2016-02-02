from oldowan.mtconvert import seq2sites, sites2seq
from string import translate
import pandas as pd
import numpy as np
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

with open('brehm2002.txt', 'rU') as f:
	data = f.readlines()

iids = []
counts = np.zeros((len(data), 9), dtype=np.int)
sites = []

for l in data:
	parts = l.strip().split()
	iids.append(parts[0])
	row = len(iids) - 1
	counts[row,] = [int(x) for x in parts[1:10]]
	# table includes some sites outside region, but are in parens
	sites_arr = [x for x in parts[10:] if not x.startswith('(')]
	sites.append(' '.join(sites_arr))

## Validate
passed_validation = True

for i in range(len(sites)):
	curr_sites = sites[i]
	seq1 = sites2seq(curr_sites, region, add16k=True)
	mysites = seq2sites(seq1)
	if not len(curr_sites.split()) == len(mysites):
		passed_validation = False
		print "different number of sites"
		print iids[i]
	else:
		mysites_positions = np.array([x.position-16000 for x in mysites])
		curr_sites_positions = np.array([int(x[:3]) for x in curr_sites.split()])
		diff = sum(mysites_positions - curr_sites_positions)
		if diff != 0:
			passed_validation = False
			print "different sites"
			print iids[i]

if passed_validation:
	counters = [1] * 9
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			curr_sites = sites[i]
			seq1 = sites2seq(curr_sites, region, add16k=True)
			mysites = seq2sites(seq1)
			mysites = ' '.join([str(x) for x in mysites])
			for j in range(9):
				if counts[i,j] > 0:
					prefix = metadata.ix[j,'NewPrefix']
					for k in range(counts[i,j]):
						newid = prefix + str(counters[j]).zfill(3)
						f.write('%s,%s,%s\n' % (newid, iids[i], mysites))
						counters[j] += 1