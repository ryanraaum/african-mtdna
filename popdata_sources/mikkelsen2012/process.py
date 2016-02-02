from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

with open('mikkelsen2012.csv', 'rU') as f:
        data = f.readlines()

hids = []
sites = []

for i in range(len(data)):
        x = data[i].strip().split(',')
        hids.append(x[0])
        sites.append(' '.join(x[1:]))

## Validate
passed_validation = True

for i in range(len(sites)):
	curr_sites = str2sites(sites[i])
	seq = sites2seq(curr_sites, region)
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	with open('processed.csv', 'w') as f:
		prefix = metadata.ix[0,'NewPrefix']
		for i in range(len(sites)):
			hid = hids[i]
			curr_sites = str2sites(sites[i])
			seq = sites2seq(curr_sites, region)
			mysites = ' '.join([str(x) for x in seq2sites(seq)])
			num = hid[3:]
			newid = prefix + num
			f.write('%s,%s,%s\n' % (newid, hid, mysites))