from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import numpy as np
import sys
import csv

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)

hids = []
sites = []

with open('alabri2012.csv', 'rU') as f:
	reader = csv.reader(f)
	reader.next() # skip past header
	for row in reader:
		hids.append(row[0])
		sites.append(str2sites(row[4]))

## Validate
passed_validation = True

for i in range(len(sites)):
	region = range2region(metadata.ix[hids[i][:2],'SeqRange'])
	seq = translate(sites2seq(sites[i], region), None, '-')
	mysites = seq2sites(seq)
	if not mysites == sites[i]:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, hids[i]

counter = {}
for k in metadata.index:
	counter[k] = 0

if passed_validation:
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			key = hid[:2]
			prefix = metadata.ix[key,'NewPrefix']
			counter[key] += 1
			num = str(counter[key]).zfill(3)
			newid = prefix + num
			seq = sites2seq(sites[i], region)
			mysites = ' '.join([str(x) for x in seq2sites(seq)])
			f.write('%s,%s,%s\n' % (newid, hid, mysites))