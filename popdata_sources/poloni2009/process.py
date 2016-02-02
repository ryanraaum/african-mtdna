from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from oldowan.fasta import fasta
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
regionstr = metadata.ix[0,'SeqRange']
regionparts = regionstr.split(';')
region1 = range2region(regionparts[0])
region2 = range2region(regionparts[1]) + range2region(regionparts[2])


ff = fasta('Poloni2009HVR1.fasta', 'r')
data1 = ff.readentries()
ff.close()

ff = fasta('Poloni2009HVR2.fasta', 'r')
data2 = ff.readentries()
ff.close()

hvr1 = {}
hvr2 = {}

for e in data1:
	k = e['name'].split()[0]
	hvr1[k] = e['sequence']

for e in data2:
	k = e['name'].split()[0]
	hvr2[k] = e['sequence']

# retaining only those with both hvr1 and hvr2
k1 = set(hvr1.keys())
k2 = set(hvr2.keys())
hids = list(k1 & k2)

sites = {}
for k in hids:
	sites[k] = seq2sites(hvr1[k])
	sites[k] = sites[k] + seq2sites(hvr2[k])

## Validate
passed_validation = True

for i in range(len(sites)):
	hid = hids[i]
	seq1 = translate(sites2seq(sites[hid], region1), None, '-')
	seq2 = translate(sites2seq(sites[hid], region2), None, '-')
	if not seq1 == hvr1[hid] and seq2 == hvr2[hid]:
		passed_validation = False
		print i, hids[i]

if passed_validation:
	counter = {}
	for k in metadata.index:
		counter[k] = 0
	with open('processed.csv', 'w') as f:
		for i in range(len(hids)):
			hid = hids[i]
			grp = hid[:1]
			mysites = ' '.join([str(x) for x in sites[hid]])
			prefix = metadata.ix[grp,'NewPrefix']
			counter[grp] += 1
			num = str(counter[grp]).zfill(3)
			newid = prefix + num
			f.write('%s,%s,%s\n' % (newid, hid, mysites))