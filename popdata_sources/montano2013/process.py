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

with open('montano2013.csv', 'rU') as f:
	f.readline() # skip past header
	data = f.readlines()

counts = np.zeros((len(data), 5), dtype=np.int)

hids = []
sites = []

for i in range(len(data)):
        x = data[i].strip().split(',')
        hids.append(x[0])
        sites.append(x[2])
        count = x[4:]
        for j in range(5):
        	if count[j] == '':
        		count[j] = '0'
        counts[i,] = [int(y) for y in count]

## Validate
passed_validation = True

for i in range(len(sites)):
	curr_sites = str2sites(sites[i])
	seq = sites2seq(curr_sites, region)
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	counter = [0] * 5
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			curr_sites = str2sites(sites[i])
			seq = sites2seq(curr_sites, region)
			mysites = ' '.join([str(x) for x in seq2sites(seq)])
			for j in range(5):
				prefix = metadata.ix[j,'NewPrefix']
				for k in range(counts[i,j]):
					counter[j] += 1
					num = str(counter[j]).zfill(3)
					newid = prefix + num
					f.write('%s,%s,%s\n' % (newid, hid, mysites))