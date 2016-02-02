from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys
import csv

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0, 'SeqRange'])

hids = []
groups = []
sites = []

with open('ely2006.csv', 'rU') as f:
	reader = csv.reader(f)
	reader.next() # skip past header
	for row in reader:
		hids.append(row[0])
		groups.append(row[1])
		sites.append(str2sites(row[3], add16k=True))

## Validate
passed_validation = True

for i in range(len(sites)):
	seq = sites2seq(sites[i], region)
	mysites = seq2sites(seq)
	if not mysites == sites[i]:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not translate(seq, None, '-') == myseq:
			passed_validation = False
			print i, hids[i]

counter = {}
for k in metadata.index:
	counter[k] = 0
	
if passed_validation:
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			key = groups[i]
			counter[key] = counter[key] + 1
			newid = metadata.ix[key,'NewPrefix'] + str(counter[key]).zfill(3)
			seq = sites2seq(sites[i], range2region(metadata.ix[key,'SeqRange']))
			seq = translate(seq, None, '-')
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			f.write('%s,%s,%s\n' % (newid, hids[i], mysites))