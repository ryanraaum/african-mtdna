from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

with open('marks2014_haplotypes.csv', 'rU') as f:
	f.readline() # skip past header
	data = f.readlines()

hids = []
sites = []

for l in data:
	parts = l.strip().split(',')
	hids.append(parts[0])
	s = ' '.join([x for x in parts[4:] if len(x) > 0])
	sites.append(s)

## Validate
passed_validation = True

for i in range(len(sites)):
	curr_sites = sites[i]
	seq = translate(sites2seq(curr_sites, region), None, '-')
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	counts = pd.read_csv('marks2014_counts.csv', index_col=0)
	counter = [0] * 6
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			curr_sites = sites[i]
			seq = translate(sites2seq(curr_sites, region), None, '-')
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			for j in range(len(metadata.index)):
				prefix = metadata.ix[metadata.index[j],'NewPrefix']
				for k in range(int(counts.ix[hid, metadata.index[j]])):
					counter[j] += 1
					num = str(counter[j]).zfill(3)
					newid = prefix + num
					f.write('%s,%s,%s\n' % (newid, hid, mysites))