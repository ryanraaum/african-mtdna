from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys
import csv

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region('1-16569')

hids = []
data = []
sites = []
counts = []
popnames = None

with open('batini2011.csv', 'rU') as f:
	reader = csv.reader(f)
	header = reader.next()
	popnames = header[3:13]
	for row in reader:
		hids.append(row[0])
		data.append(row[13].split(','))
		counts.append(row[3:13])

# convert counts to integers
newcounts = []
for i in range(len(counts)):
	newcounts.append([int(x) for x in counts[i]])

for i in range(len(data)):
	sites.append(str2sites(' '.join(data[i])))

## Validate variant sites
passed_validation = True

for i in range(len(sites)):
	seq = sites2seq(sites[i], region).upper()
	mysites = seq2sites(seq)
	if not sites[i] == mysites:
		if not translate(seq, None, '-') == translate(sites2seq(mysites, region), None, '-'):
			passed_validation = False
			print i

counter = {}
for k in popnames:
	counter[k] = 0

if passed_validation:
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			seq = sites2seq(sites[i], region)
			seq = translate(seq, None, '-')
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			for j in range(len(popnames)):
				for k in range(newcounts[i][j]):
					counter[popnames[j]] += 1
					newid = metadata.ix[popnames[j], 'NewPrefix'] + str(counter[popnames[j]]).zfill(3)
					f.write('%s,%s,%s\n' % (newid, hids[i], mysites))