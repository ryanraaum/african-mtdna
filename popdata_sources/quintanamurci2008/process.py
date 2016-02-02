from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

with open('quintanamurci2008.csv', 'rU') as f:
	f.readline() # skip past header
	data = f.readlines()

grps = []
hids = []
sites = []

for l in data:
	e = l.strip().split(',')
	grps.append(e[1])
	hids.append(e[4])
	sites.append(e[5])

## Validate
passed_validation = True

for i in range(len(sites)):
	curr_sites = str2sites(sites[i])
	seq = translate(sites2seq(curr_sites, region), None, '-').upper()
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i


if passed_validation:
	counter = {}
	for k in metadata.index:
		counter[k] = 0
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			grp = grps[i]
			hid = hids[i]
			prefix = metadata.ix[grp,'NewPrefix']
			seq = translate(sites2seq(sites[i], region), None, '-') 
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			counter[grp] += 1
			num = str(counter[grp]).zfill(3)
			newid = prefix + num
			f.write('%s,%s,%s\n' % (newid, hid, mysites))