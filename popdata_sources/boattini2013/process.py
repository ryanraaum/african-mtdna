from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import re
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

with open('boattini2013.csv', 'rU') as f:
	data = f.readlines()

## Validate variant sites
passed_validation = True

for l in data:
	parts = l.strip().split(',')
	sites = parts[1].split()
	sites.sort()
	sites = ' '.join(sites)
	seq1 = sites2seq(parts[1], region)
	mysites = seq2sites(seq1)
	if not sites == ' '.join([str(x) for x in mysites]):
		if not translate(seq1, None, '-') == translate(sites2seq(mysites, region), None, '-'):
			passed_validation = False
			print l

if passed_validation:
	with open('processed.csv', 'w') as f:
		for l in data:
			parts = l.strip().split(',')
			origid = parts[0]
			key = parts[0][:3]
			m = re.search(r'[0-9]+', parts[0])
			counter = m.group()
			newid = metadata.ix[key,'NewPrefix'] + counter.zfill(3)
			sites = parts[1].split()
			seq = sites2seq(sites, range2region(metadata.ix[key,'SeqRange']))
			seq = translate(seq, None, '-')
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			f.write('%s,%s,%s\n' % (newid, origid, mysites))