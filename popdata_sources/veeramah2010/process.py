from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import numpy as np
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)

with open('veeramah2010.csv', 'rU') as f:
	header = f.readline() 
	data = f.readlines()

header = header.strip().split(',')
popnames = [x.strip() for x in header[2:]]
counts = np.zeros((len(data), 34), dtype=np.int)

hids = []
sites = []

for i in range(len(data)):
        x = data[i].strip().split('"')
        hids.append(x[0].split(',')[0])
        y = x[1].strip().split(',')
        y = [a.strip() for a in y if len(a) > 0]
        sites.append(' '.join(y))
        count = x[2].split(',')[1:]
        for j in range(counts.shape[1]):
        	if count[j] == '':
        		count[j] = '0'
        counts[i,] = [int(y) for y in count]

counts = pd.DataFrame(counts, columns=popnames)

## Validate
passed_validation = True

# use larger region for validation
region = range2region('16000-16400')

for i in range(len(sites)):
	x = sites[i]
	curr_sites = str2sites(x, add16k=True)
	seq = sites2seq(curr_sites, region)
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	counter = pd.Series([0] * counts.shape[1], index=counts.columns)
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			for pop in metadata.index:
				prefix = metadata.ix[pop,'NewPrefix']
				region = range2region('16000-16400')
				s = sites[i]
				curr_sites = str2sites(s, add16k=True)
				seq = sites2seq(curr_sites, region)
				mysites = ' '.join([str(x) for x in seq2sites(seq)])
				for k in range(counts.ix[i,pop]):
					counter[pop] += 1
					num = str(counter[pop]).zfill(3)
					newid = prefix + num
					f.write('%s,%s,%s\n' % (newid, hid, mysites))