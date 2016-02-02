from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from oldowan.polymorphism import Polymorphism
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

with open('gonzalez2006_haplotypes.csv', 'rU') as f:
	f.readline() # skip past header
	data = f.readlines()

hids = []
sites = []

for l in data:
	parts = l.strip().split(',')
	hids.append(parts[0])
	sites.append('%s %s' % (parts[1],parts[2]))

## need to preprocess sites data because Gonzalez et al. use some nonstandard notation
newsites = []
for i in range(len(sites)):
	s = str2sites(sites[i])
	newsites.append([])
	for j in range(len(s)):
		val = s[j].value
		if val.startswith('d'):
			s[j].value = '-'
			newsites[i].append(s[j])
			if len(val) > 2:
				p = Polymorphism(s[j].position+1,0,'-')
				newsites[i].append(p)
		elif val.startswith('i'):
			if val.startswith('ii'):
				val = val[1:]
			inserts = list(val[1:])
			for k in range(len(inserts)):
				p = Polymorphism(s[j].position, k+1, inserts[k])
				newsites[i].append(p)
		elif val.startswith('.'):
			pos = s[j].position
			p = Polymorphism(pos, 1, sites2seq('', (pos,pos)))
			newsites[i].append(p)
		else:
			newsites[i].append(s[j])

## Validate
passed_validation = True

for i in range(len(newsites)):
	curr_sites = newsites[i]
	seq = translate(sites2seq(curr_sites, region), None, '-')
	mysites = seq2sites(seq)
	if not mysites == curr_sites:
		myseq = translate(sites2seq(mysites, region), None, '-')
		if not seq == myseq:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	counts = pd.read_csv('gonzalez2006_counts.csv', index_col=0)
	counter = [1] * 9
	with open('processed.csv', 'w') as f:
		for i in range(len(newsites)):
			hid = hids[i]
			seq = translate(sites2seq(newsites[i], region), None, '-')
			mysites = ' '.join([str(x) for x in seq2sites(seq)])
			for j in range(len(metadata.index)):
				prefix = metadata.ix[metadata.index[j],'NewPrefix']
				for k in range(counts.ix[hid, metadata.index[j]]):
					num = str(counter[j]).zfill(3)
					counter[j] += 1
					newid = prefix + num
					f.write('%s,%s,%s\n' % (newid, hid, mysites))