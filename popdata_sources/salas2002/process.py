from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

counts = pd.read_table('salas2002_counts.txt', index_col=0, sep='\W')

with open('salas2002_haplotypes.txt', 'rU') as f:
	data = f.readlines()

hids = []
sites = []

for l in data:
	e = l.strip().split()
	hids.append(e[0])
	sites.append(' '.join(e[1:]))

## Validate
passed_validation = True

for i in range(len(sites)):
	curr_sites = str2sites(sites[i], add16k=True)
	seq = translate(sites2seq(curr_sites, region), None, '-')
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	counter = {}
	for k in metadata.index:
		counter[k] = 0
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			seq = translate(sites2seq(sites[i], region, add16k=True), None, '-') 
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			for j in metadata.index:
				for k in range(counts.ix[hid,j]):
					prefix = metadata.ix[j,'NewPrefix']
					counter[j] += 1
					num = str(counter[j]).zfill(3)
					newid = prefix + num
					f.write('%s,%s,%s\n' % (newid, hid, mysites))