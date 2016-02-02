from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import re
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
regionstr = metadata.ix[0,'SeqRange']
regionparts = regionstr.split(';')
region1 = range2region(regionparts[0])
region2 = range2region(regionparts[1])

with open('pereira2001.txt', 'rU') as f:
	data = f.readlines()

freq = []
hvr1 = []
hvr2 = []
hgrp = []

for l in data:
	e = l.strip().split(',')
	freq.append(int(e[0]))
	hvr1.append(e[1])
	hvr2.append(e[2])
	hgrp.append(e[3])

# some sites have nonstandard format, fix those
for i in range(len(freq)):
	h1 = hvr1[i].split()
	h1f = []
	for s in h1:
		if '/' in s:
			m = re.match(r'([0-9]+)[A-Z]/([A-Z])', s)
			s = '%s%s' % m.groups()
		h1f.append(s)
	hvr1[i] = ' '.join(h1f)

	h2 = hvr2[i].split()
	h2f = []
	for s in h2:
		if '/' in s:
			m = re.match(r'([0-9]+)[A-Z]/([A-Z])', s)
			s = '%s%s' % m.groups()
		elif '.1' in s:
			s = s + 'C'
		elif '.2' in s:
			h2f.append(s[:-1] + '1C')
			s = s + 'C'
		h2f.append(s)
	hvr2[i] = ' '.join(h2f)

## Validate
passed_validation = True

for i in range(len(freq)):
	curr_sites = str2sites(hvr1[i], add16k=True)
	seq = translate(sites2seq(curr_sites, region1), None, '-')
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region1), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, 'hvr1'
	curr_sites = str2sites(hvr2[i])
	seq = translate(sites2seq(curr_sites, region2), None, '-')
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region2), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, 'hvr2'


if passed_validation:
	counter = 0
	prefix = metadata.ix[0,'NewPrefix']
	with open('processed.csv', 'w') as f:
		for i in range(len(freq)):
			hid = hgrp[i]
			seq1 = translate(sites2seq(hvr1[i], region1, add16k=True), None, '-') 
			seq2 = translate(sites2seq(hvr2[i], region2), None, '-') 
			mysites = seq2sites(seq1) + seq2sites(seq2)
			mysites = ' '.join([str(x) for x in mysites])
			for j in range(freq[i]):
				counter += 1
				num = str(counter).zfill(3)
				newid = prefix + num
				f.write('%s,%s,%s\n' % (newid, hid, mysites))