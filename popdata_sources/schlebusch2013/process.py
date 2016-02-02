from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import re
import sys

sys.path.append('../../scripts')
from utils import *
from genbank import read_genbank

def getnote(e):
	for x in e['features']:
		if x[0] == 'source':
			for y in x[1]:
				if isinstance(y, tuple):
					if y[0] == 'note':
						return y[1]
	return None

def getisolate(e):
	for x in e['features']:
		if x[0] == 'source':
			for y in x[1]:
				if isinstance(y, tuple):
					if y[0] == 'isolate':
						return y[1]
	return None


## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
regionparts = metadata.ix[0, 'SeqRange'].split(';')
region1 = range2region(regionparts[0])
region2 = range2region(regionparts[1])

entries1 = read_genbank('hvr1.gb', what='filename')

hids1 = []
groups1 = []
seqs1 = []
sites1 = []

for e in entries1:
	hid = getisolate(e)
	if hid is not None:
		hids1.append(hid)
	n = getnote(e)
	if n is not None:
		parts = n.split()
		groups1.append(parts[-1])
	seqs1.append(e['sequence'])
	sites1.append(seq2sites(e['sequence']))

entries2 = read_genbank('hvr2.gb', what='filename')

hids2 = []
groups2 = []
seqs2 = []
sites2 = []

for e in entries2:
	hid = getisolate(e)
	if hid is not None:
		hids2.append(hid)
	n = getnote(e)
	if n is not None:
		parts = n.split()
		groups2.append(parts[-1])
	seqs2.append(e['sequence'])
	sites2.append(seq2sites(e['sequence']))

prefixes = []
for h in hids1:
	m = re.match(r'([a-zA-Z]+)[0-9]*', h)
	if m is not None:
		prefixes.append(m.groups()[0])
	else:
		print h

## Validate
passed_validation = True

if hids1 == hids2:
	hids = hids1
else:
	passed_validation = False
	print "Haplotype ids don't match across HVR1 and HVR2"

if groups1 == groups2:
	groups = groups1
else:
	passed_validation = False
	print "Group ids don't match across HVR1 and HVR2"

#combine sites
sites = []
for i in range(len(hids)):
	sites.append(sites1[i] + sites2[i])

for i in range(len(hids)):
	hid = hids[i]
	seq1 = translate(sites2seq(sites[i], region1), None, '-')
	seq2 = translate(sites2seq(sites[i], region2), None, '-')
	if not seq1 == seqs1[i].upper() and seq2 == seqs2[i].upper():
		passed_validation = False
		print i, hids[i]

if passed_validation:
	counter = {}
	for k in metadata.index:
		counter[k] = 0
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			grp = groups[i]
			mysites = ' '.join([str(x) for x in sites[i]])
			prefix = metadata.ix[grp,'NewPrefix']
			counter[grp] += 1
			num = str(counter[grp]).zfill(3)
			newid = prefix + num
			f.write('%s,%s,%s\n' % (newid, hid, mysites))