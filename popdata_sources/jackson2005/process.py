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

with open('jackson2005.txt', 'rU') as f:
	data = f.readlines()

hgrp = []
counts = np.zeros((len(data), 4), dtype=np.int)
sites = []

for l in data:
	e = l.strip().split()
	counts[len(hgrp),] = [int(x) for x in e[1:5]]
	hgrp.append(e[0])
	sites.append(' '.join(e[5:]))

## Validate
passed_validation = True

for i in range(len(sites)):
	curr_sites = str2sites(sites[i], add16k=True)
	seq = translate(sites2seq(curr_sites, region), None, '-')
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	counter = [0] * 4
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hgrp[i]
			seq = translate(sites2seq(sites[i], region, add16k=True), None, '-')
			mysites = ' '.join([str(x) for x in seq2sites(seq)])
			for j in range(len(metadata.index)):
				prefix = metadata.ix[metadata.index[j],'NewPrefix']
				for k in range(counts[i, j]):
					counter[j] += 1
					num = str(counter[j]).zfill(3)
					newid = prefix + num
					f.write('%s,%s,%s\n' % (newid, hid, mysites))