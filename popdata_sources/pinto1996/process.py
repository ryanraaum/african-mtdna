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

with open('pinto1996.csv', 'rU') as f:
        header = f.readline()
        reference = f.readline()
        data = f.readlines()

positions = header.strip().split(',')[1:-4]
reference = reference.strip().split(',')[1:-4]

counts = np.zeros((len(data), 3), dtype=np.int)

hids = []
vals = []

for i in range(len(data)):
        x = data[i].strip().split(',')
        hids.append(x[0])
        vals.append(x[1:-4])
        counts[i,] = [int(x) for x in x[-3:]]

# numbering is weird. 23 = 16019, so need to subtract 4 from header values
positions = [int(x) - 4 + 16000 for x in positions]
positions = [str(x) for x in positions]

sites = []
for i in range(len(vals)):
	x = vals[i]
	y = []
	for j in range(len(x)):
		if x[j] != '.':
			y.append('%s%s' % (positions[j], x[j]))
	sites.append(' '.join(y))

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
	counter = [0] * 3 
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			curr_sites = str2sites(sites[i])
			seq = sites2seq(curr_sites, region)
			mysites = ' '.join([str(x) for x in seq2sites(seq)])
			for j in range(3):
				prefix = metadata.ix[j,'NewPrefix']
				for k in range(counts[i,j]):
					counter[j] += 1
					num = str(counter[j]).zfill(3)
					newid = prefix + num
					f.write('%s,%s,%s\n' % (newid, hid, mysites))