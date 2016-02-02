from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)

regionparts = metadata.ix[0,'SeqRange'].split(';')
region1 = range2region(regionparts[0])
region2 = range2region(regionparts[1])

with open('coelho2009_haplotypes.csv', 'rU') as f:
	f.readline() # skip past header
	data = f.readlines()

hids = []
hvr1sites = []
hvr2sites = []

for l in data:
	parts = l.strip().split(',')
	if int(parts[3]) == 377 and int(parts[7]) == 268: 
		hids.append(parts[0])
		hvr1sites.append(parts[4])
		hvr2sites.append(parts[8])

## need to preprocess sites data for some nonstandard notation in hvr2
hvr1 = []
hvr2 = []
for i in range(len(hids)):
	s1 = str2sites(hvr1sites[i], add16k=True)
	hvr1.append(s1)

	s2 = hvr2sites[i].split()
	s2new = []
	for j in range(len(s2)):
		if s2[j].endswith('.2C'):
			parts = s2[j].split('.')
			s2new.append('%s.1C' % parts[0])
			s2new.append('%s.2C' % parts[0])
		else:
			s2new.append(s2[j])
	s2 = str2sites(' '.join(s2new))
	hvr2.append(s2)

newsites = []
for i in range(len(hvr1)):
	newsites.append(hvr1[i] + hvr2[i])

## Validate
passed_validation = True

for i in range(len(newsites)):
	curr_sites = newsites[i]
	seq1 = translate(sites2seq(curr_sites, region1), None, '-')
	seq2 = translate(sites2seq(curr_sites, region2), None, '-')
	mysites = seq2sites(seq1) + seq2sites(seq1)
	if not mysites == curr_sites:
		myseq1 = translate(sites2seq(mysites, region1), None, '-')
		myseq2 = translate(sites2seq(mysites, region2), None, '-')
		if not seq1 == myseq1 and seq2 == myseq2:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	counts = pd.read_csv('coelho2009_counts.csv', index_col=0)
	counts = counts.fillna(0)
	counter = [0] * 5
	with open('processed.csv', 'w') as f:
		for i in range(len(newsites)):
			hid = hids[i]
			curr_sites = newsites[i]
			seq1 = translate(sites2seq(curr_sites, region1), None, '-')
			seq2 = translate(sites2seq(curr_sites, region2), None, '-')
			mysites = seq2sites(seq1) + seq2sites(seq2)
			mysites = ' '.join([str(x) for x in mysites])
			for j in range(len(metadata.index)):
				prefix = metadata.ix[metadata.index[j],'NewPrefix']
				for k in range(int(counts.ix[hid, metadata.index[j]])):
					counter[j] += 1
					num = str(counter[j]).zfill(3)
					newid = prefix + num
					f.write('%s,%s,%s\n' % (newid, hid, mysites))