from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

with open('kivisild2004.csv', 'rU') as f:
	f.readline() # skip past header (3 lines)
	f.readline() 
	f.readline()
	data = f.readlines()

accn = []
group = []
sites = []

for l in data:
	e = l.strip().split(',')
	accn.append(e[1])
	group.append(e[2])
	sites.append(e[4])

# some sites are outside the specified region and are in parens, strip those out
newsites = []
for s in sites:
	x = s.split()
	x = [y for y in x if not y.startswith('(')]
	newsites.append(' '.join(x))

## Validate
passed_validation = True

for i in range(len(newsites)):
	curr_sites = str2sites(newsites[i], add16k=True)
	seq = translate(sites2seq(curr_sites, region), None, '-')
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, accn[i]

if passed_validation:
	counter = {}
	for k in metadata.index:
		counter[k] = 0
	with open('processed.csv', 'w') as f:
		for i in range(len(newsites)):
			hid = accn[i]
			grp = group[i]
			seq = translate(sites2seq(newsites[i], region, add16k=True), None, '-')
			mysites = ' '.join([str(x) for x in seq2sites(seq)])
			if grp == 'ERI':
				grp = 'TIG' # the Eritrean labeled individuals are Tigray
			if grp not in metadata.index:
				grp = 'OTH'
			prefix = metadata.ix[grp,'NewPrefix']
			counter[grp] += 1
			num = str(counter[grp]).zfill(3)
			newid = prefix + num
			f.write('%s,%s,%s\n' % (newid, hid, mysites))