from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import numpy as np
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

with open('watson1997.csv', 'rU') as f:
	header = f.readline()
	data = f.readlines()

popnames = header.strip().split(',')[1:]
counts = pd.DataFrame(np.zeros((len(data), 13), dtype=np.int), columns=popnames)
sites = []

for i in range(len(data)):
	parts = data[i].strip().split(',')
	sites.append(parts[0])
	count = [int(x) for x in parts[1:]]
	counts.ix[i,] = count

# some non-canonical site names to be fixed
newsites = []

for i in range(len(sites)):
	s = sites[i].split()
	s2 = []
	for x in s:
		if 'G/A' in x:
			x = x[:3] + 'R'
		elif '-' in x:
			x = x[:3] + x[-1]
		s2.append(x)
	newsites.append(' '.join(s2))

## Validate
passed_validation = True

for i in range(len(newsites)):
	curr_sites = str2sites(newsites[i], add16k=True)
	seq = sites2seq(curr_sites, region)
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i

if passed_validation:
	counter = pd.Series([0] * counts.shape[1], index=counts.columns)
	with open('processed.csv', 'w') as f:
		for i in range(len(newsites)):
			s = newsites[i]
			curr_sites = str2sites(s, add16k=True)
			seq = sites2seq(curr_sites, region)
			mysites = ' '.join([str(x) for x in seq2sites(seq)])
			for pop in metadata.index:
				prefix = metadata.ix[pop,'NewPrefix']
				for k in range(counts.ix[i,pop]):
					counter[pop] += 1
					num = str(counter[pop]).zfill(3)
					newid = prefix + num
					f.write('%s,%s,%s\n' % (newid, str(i+1), mysites))