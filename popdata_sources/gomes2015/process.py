from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys
import csv

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)

groups = []
hids = []
hvr1 = []
hvr2 = []
sites = []

with open('gomes2015.csv', 'rU') as f:
	reader = csv.reader(f)
	reader.next() # skip past header
	for row in reader:
		groups.append(row[3])
		hids.append(row[1])
		hvr1.append(str2sites(row[6], add16k=True))
		hvr2.append(str2sites(row[8]))

for i in range(len(groups)):
	sites.append(hvr1[i] + hvr2[i])

## Validate variant sites
passed_validation = True

for i in range(len(groups)):
	region = range2region(metadata.ix[groups[i],'SeqRange'])
	seq = sites2seq(sites[i], region)
	mysites = seq2sites(seq)
	if not sites[i] == mysites:
		if not translate(seq, None, '-') == translate(sites2seq(mysites, region), None, '-'):
			passed_validation = False
			print i

counter = {}
for k in metadata.index:
	counter[k] = 0

if passed_validation:
	with open('processed.csv', 'w') as f:
		for i in range(len(groups)):
			key = groups[i]
			counter[key] = counter[key] + 1
			newid = metadata.ix[key,'NewPrefix'] + str(counter[key]).zfill(3)
			seq = sites2seq(sites[i], range2region(metadata.ix[key,'SeqRange']))
			seq = translate(seq, None, '-')
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			f.write('%s,%s,%s\n' % (newid, hids[i], mysites))