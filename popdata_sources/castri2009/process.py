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
sh = []
rw = []

with open('castri2009.csv', 'rU') as f:
	reader = csv.reader(f)
	reader.next() # skip past header
	for row in reader:
		hids.append(row[0])
		sites.append(str2sites(row[1], add16k=True))
		sh.append(int(row[3]))
		rw.append(int(row[4]))

## Validate
passed_validation = True

for i in range(len(sites)):
	seq = sites2seq(sites[i], region)
	mysites = seq2sites(seq)
	if not mysites == sites[i]:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	shcount = 0
	rwcount = 0
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			seq = sites2seq(sites[i], region)
			mysites = ' '.join([str(x) for x in seq2sites(seq)])
			for j in range(sh[i]):
				prefix = metadata.ix['Shona','NewPrefix']
				shcount = shcount + 1
				newid = prefix + str(shcount).zfill(3)
				f.write('%s,%s,%s\n' % (newid, hids[i], mysites))
			for j in range(rw[i]):
				prefix = metadata.ix['Hutu','NewPrefix']
				rwcount = rwcount + 1
				newid = prefix + str(rwcount).zfill(3)
				f.write('%s,%s,%s\n' % (newid, hids[i], mysites))