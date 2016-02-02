from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys
import csv

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

counts = pd.read_csv('plaster2011_counts.csv', index_col=0)
counts = counts.fillna(0)

hids = []
sites = []

with open('plaster2011_haplotypes.csv', 'rU') as f:
	reader = csv.reader(f)
	reader.next() # skip past header
	for row in reader:
		hids.append(row[0])
		parts = row[1].split(',')
		sites.append(str2sites(' '.join(parts), add16k=True))

## Validate
passed_validation = True

for i in range(len(hids)):
	seq = sites2seq(sites[i], region)
	mysites = seq2sites(seq)
	if not sites[i] == mysites:
		if not translate(seq, None, '-') == translate(sites2seq(mysites, region), None, '-'):
			passed_validation = False
			print i, hids[i]

counter = {}
for k in counts.columns:
	counter[k] = 0

if passed_validation:
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			seq = translate(sites2seq(sites[i], region), None, '-')
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			for key in counts.columns:
				prefix = metadata.ix[key,'NewPrefix']
				for j in range(int(counts.ix[hid, key])):
					counter[key] += 1
					num = str(counter[key]).zfill(3)
					newid = prefix + num
					f.write('%s,%s,%s\n' % (newid, hid, mysites))