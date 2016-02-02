from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

region_str = metadata.ix[0,'SeqRange']
region_parts = region_str.split(';')
region1 = range2region(region_parts[0])
region2 = range2region(region_parts[1])

with open('graven1995.csv', 'rU') as f:
	header = f.readline()
	reference = f.readline()
	data = f.readlines()

positions = header.strip().split(',')[2:]
reference = reference.strip().split(',')[2:]

hids = []
freq = []
vals = []

for l in data:
	x = l.strip().split(',')
	hids.append(x[0])
	freq.append(int(x[1]))
	vals.append(x[2:])

sites = []
for i in range(len(vals)):
	x = vals[i]
	y = []
	for j in range(len(x)):
		if x[j] != '.':
			if j == 83:
				y.append('309.1C')
			elif j == 84:
				y.append('309.2C')
			elif j == 85:
				y.append('313.1C')
			elif j == 86:
				y.append('315.1C')
			else:
				y.append('%s%s' % (positions[j], x[j]))
	sites.append(' '.join(y))

## Validate
passed_validation = True

for i in range(len(sites)):
	curr_sites = str2sites(sites[i])
	cseq1 = sites2seq(curr_sites, region1)
	cseq2 = sites2seq(curr_sites, region2)
	mysites1 = seq2sites(cseq1)
	mysites2 = seq2sites(cseq2)
	mysites = mysites1 + mysites2
	if not mysites == curr_sites:
		seq = cseq1 + cseq2
		myseq = translate(sites2seq(mysites, region1), None, '-') + translate(sites2seq(mysites, region2), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	count = 0
	prefix = metadata.ix[0,'NewPrefix']
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			curr_sites = str2sites(sites[i])
			cseq1 = sites2seq(curr_sites, region1)
			cseq2 = sites2seq(curr_sites, region2)
			mysites1 = seq2sites(cseq1)
			mysites2 = seq2sites(cseq2)
			mysites = mysites1 + mysites2
			mysites = ' '.join([str(x) for x in mysites])
			for j in range(freq[i]):
				count += 1
				num = str(count).zfill(3)
				newid = prefix + num
				f.write('%s,%s,%s\n' % (newid, hid, mysites))