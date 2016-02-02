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
region = range2region(metadata.ix[0,'SeqRange'])

counts = []
sites = []
popnames = None

with open('rando1998.csv', 'rU') as f:
	reader = csv.reader(f)
	header = reader.next()
	popnames = header[1:]
	for row in reader:
		sites.append(str2sites(row[0], add16k=True))
		counts.append(row[1:])

countm = np.zeros((len(counts), len(popnames)), dtype=np.int)
for i in range(len(counts)):
	countm[i] = [int(x) for x in counts[i]]

## Validate
passed_validation = True

for i in range(len(sites)):
	seq = sites2seq(sites[i], region)
	mysites = seq2sites(seq)
	if not mysites == sites[i]:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not translate(seq, None, '-') == myseq:
			passed_validation = False
			print i

counter = {}
for name in popnames:
	counter[name] = 0

if passed_validation:
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = "H%d" % (i + 1,)
			seq = sites2seq(sites[i], region)
			mysites = ' '.join([str(x) for x in seq2sites(seq)])
			for j in range(len(popnames)):
				prefix = metadata.ix[j,'NewPrefix']
				for k in range(countm[i,j]):
					counter[popnames[j]] += 1
					num = str(counter[popnames[j]]).zfill(3)
					newid = prefix + num
					f.write('%s,%s,%s\n' % (newid, hid, mysites))