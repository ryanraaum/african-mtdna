from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys
import csv

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region('16024-16569;1-400')

hids = []
hvr1 = []
hvr2 = []
sites = []

with open('soares2012.csv', 'rU') as f:
	reader = csv.reader(f)
	reader.next() # skip past header
	for row in reader:
		for i in range(int(row[1])):
			hids.append(row[0])
			hvr1.append(row[2].strip().split())
			hvr2.append(row[3].strip().split())

# transversions are reported as A/D (rcrs/derived), which is not standard notation
# so need to go through and fix all of those
def fix(x):
	if '/' in x:
		derived = x[-1]
		x = x[:-3] + x[-1]
	return x	

for i in range(len(hvr1)):
	hvr1[i] = [fix(x) for x in hvr1[i]]

for i in range(len(hvr2)):
	hvr2[i] = [fix(x) for x in hvr2[i]]

for i in range(len(hids)):
	sites.append(str2sites(' '.join(hvr1[i]), add16k=True) + str2sites(' '.join(hvr2[i])))

## Validate variant sites
passed_validation = True

for i in range(len(hids)):
	seq = sites2seq(sites[i], region)
	mysites = seq2sites(seq)
	if not sites[i] == mysites:
		if not translate(seq, None, '-') == translate(sites2seq(mysites, region), None, '-'):
			passed_validation = False
			print i, hids[i]

counter = {}
for k in metadata.index:
	counter[k] = 0

if passed_validation:
	with open('processed.csv', 'w') as f:
		for i in range(len(hids)):
			hid = hids[i]
			key = hid[:3]
			counter[key] += 1
			newid = metadata.ix[key,'NewPrefix'] + str(counter[key]).zfill(3)
			seq = sites2seq(sites[i], region)
			seq = translate(seq, None, '-')
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			f.write('%s,%s,%s\n' % (newid, hid, mysites))