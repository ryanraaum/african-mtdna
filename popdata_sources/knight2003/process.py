from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from oldowan.fasta import fasta
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

ff = fasta('knight2003.fasta', 'r')
data = ff.readentries()
ff.close()

hids = []
seqs = []
sites = []

for e in data:
	hids.append(e['name'].split()[0])
	seqs.append(e['sequence'])
	m = re.search(r'GATCACA', e['sequence'])
	cut = m.start()
	seq1 = e['sequence'][:cut]
	seq2 = e['sequence'][cut:]
	sites.append(seq2sites(seq1) + seq2sites(seq2))

## Validate
passed_validation = True

for i in range(len(sites)):
	seq = translate(sites2seq(sites[i], region1), None, '-') + translate(sites2seq(sites[i], region2), None, '-')
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