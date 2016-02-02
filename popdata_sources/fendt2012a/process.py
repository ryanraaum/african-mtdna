from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

with open('fendt2012.csv', 'rU') as f:
	data = f.readlines()

hids = []
sites = []

for l in data:
	parts = l.strip().split(',')
	hids.append(parts[0])
	nonempty = [x for x in parts[3:] if len(x) > 0]
	sites.append(' '.join(nonempty))

## Validate
passed_validation = True

for i in range(len(sites)):
	curr_sites = str2sites(sites[i])
	seq = translate(sites2seq(curr_sites, region), None, '-')
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	prefix = metadata.ix[0,'NewPrefix']
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			seq = translate(sites2seq(sites[i], region), None, '-')
			mysites = ' '.join([str(x) for x in seq2sites(seq)])
			num = hids[i][-3:]
			newid = prefix + num
			f.write('%s,%s,%s\n' % (newid, hids[i], mysites))