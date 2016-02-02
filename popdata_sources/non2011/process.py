from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from oldowan.fasta import fasta
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)

ff = fasta('non_2011.fasta', 'r')
data = ff.readentries()
ff.close()

hids = []
seqs = []
sites = []

# four sequences are shorter than all the rest, will drop them

for e in data:
	if len(e['sequence']) > 350:
		hids.append(e['name'].split()[0])
		seqs.append(e['sequence'])
		sites.append(seq2sites(e['sequence']))

## Validate
passed_validation = True

for i in range(len(sites)):
	hid = hids[i]
	key = hid[:2]
	region = range2region(metadata.ix[key, 'SeqRange'])
	seq = translate(sites2seq(sites[i], region), None, '-')
	if not seq == seqs[i]:
		passed_validation = False
		print i, hids[i]

if passed_validation:
	counter = {}
	for k in metadata.index:
		counter[k] = 0
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			grp = hid[:2]
			mysites = ' '.join([str(x) for x in sites[i]])
			prefix = metadata.ix[grp,'NewPrefix']
			counter[grp] += 1
			num = str(counter[grp]).zfill(3)
			newid = prefix + num
			f.write('%s,%s,%s\n' % (newid, hid, mysites))