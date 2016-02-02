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

hids = []
sites = []

with open('poetsch2013.csv', 'rU') as f:
	reader = csv.reader(f)
	reader.next() # skip past header
	for row in reader:
		hids.append(row[0])
		sitestr = ' '.join(row[1:])
		sites.append(str2sites(sitestr))


## Validate variant sites
passed_validation = True

for i in range(len(hids)):
	seq = sites2seq(sites[i], region)
	mysites = seq2sites(seq)
	if not sites[i] == mysites:
		if not translate(seq, None, '-') == translate(sites2seq(mysites, region), None, '-'):
			passed_validation = False
			print i


if passed_validation:
	counter = 0
	prefix = metadata.ix[0,'NewPrefix']
	with open('processed.csv', 'w') as f:
		for i in range(len(hids)):
			counter = counter + 1
			newid = prefix + str(counter).zfill(3)
			seq = sites2seq(sites[i], region)
			seq = translate(seq, None, '-')
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			f.write('%s,%s,%s\n' % (newid, hids[i], mysites))