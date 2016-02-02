from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys
import csv

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region('16024-16569;1-340')

hids = []
hvr1 = []
hvr2 = []
sites = []

with open('beleza2005.csv', 'rU') as f:
	reader = csv.reader(f)
	reader.next() # skip past header
	for row in reader:
		for i in range(int(row[1])):
			hids.append(row[0])
			hvr1.append(str2sites(row[2], add16k=True))
			hvr2.append(str2sites(row[3]))

for i in range(len(hids)):
	sites.append(hvr1[i] + hvr2[i])

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