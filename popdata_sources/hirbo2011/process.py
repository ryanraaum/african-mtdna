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

with open('hirbo2011.csv', 'rU') as f:
	data = f.readlines()

hids = []
sites = []

for l in data:
	parts = l.strip().split('"')
	parts = [x for x in parts if len(x) > 0]
	s = parts[1]
	s = [x.strip() for x in s.split(',')]
	h = parts[-1]
	h = [x.strip() for x in h.split(',')]
	h = [x for x in h if len(x) > 3]
	for i in h:
		hids.append(i)
		sites.append(s)

newsites = []
# fix some non-canonical sites
for row in sites:
	newrow = []
	for s in row:
		if ':' in s:
			# these are deletions
			newrow.append(s[:-1] + 'd')
		elif s == '522-3D' or s == '522- 3D':
			newrow.append('522d')
			newrow.append('523d')
		elif s == '16318CT':
			newrow.append('16318Y')
		elif '+' in s:
			parts = s.split('+')
			loc = parts[0]
			count = parts[1][0]
			nuc = parts[1][1:]
			if len(nuc) == 1:
				for i in range(int(count)):
					newrow.append('%s.%s%s' % (loc, str(i+1), nuc))
			else:
				count = int(count) * len(nuc)
				for i in range(int(count)):
					newrow.append('%s.%s%s' % (loc, str(i+1), nuc[i % len(nuc)]))
		elif '.' in s:
			parts = s.split('.')
			if re.match(r'[0-9][A-Z][A-Z]', parts[-1]) is not None:
				count = parts[1][0]
				nuc = parts[1][1:]
				count = int(count) * len(nuc)
				for i in range(int(count)):
					newrow.append('%s.%s%s' % (loc, str(i+1), nuc[i % len(nuc)]))
			else:
				newrow.append(s)
		else:
			newrow.append(s)
	newsites.append(' '.join(newrow))

## Validate
passed_validation = True

for i in range(len(newsites)):
	curr_sites = str2sites(newsites[i])
	# some entries have data outside explicitly sequenced 15900-640 region
	# get rid of extra sites
	curr_sites = [x for x in curr_sites if x.position >= 15900 or x.position <= 640]
	seq = sites2seq(curr_sites, region)
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = sites2seq(mysites, region)
		if not seq == myseq:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	counter = {}
	for i in metadata.index:
		counter[i] = 0
	with open('processed.csv', 'w') as f:
		for i in range(len(newsites)):
			curr_sites = str2sites(newsites[i])
			# some entries have data outside explicitly sequenced 15900-640 region
			# get rid of extra sites
			curr_sites = [x for x in curr_sites if x.position >= 15900 or x.position <= 640]
			seq = sites2seq(curr_sites, region)
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			for key in metadata.index:
				if hids[i].upper().startswith(key):
					prefix = metadata.ix[key, 'NewPrefix']		
					counter[key] += 1
					newid = prefix + str(counter[key]).zfill(3)						
					f.write('%s,%s,%s\n' % (newid, hids[i], mysites))
