from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from oldowan.fasta import fasta
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0, 'SeqRange'])

ff = fasta('coia2005.fasta', 'r')
data = ff.readentries()
ff.close()

hids = []
sites = []

for entry in data:
	words = entry['name'].split()
	hids.append(words[4])
	sites.append(seq2sites(entry['sequence']))

# validate
passed_validation = True

for i in range(len(sites)):
	seq1 = data[i]['sequence']
	if not seq1 == translate(sites2seq(sites[i], region), None, '-'):
			passed_validation = False
			print i, hids[i]

counter = {}
for k in metadata.index:
	counter[k] = 0

if passed_validation:
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			key = None
			if hid.startswith('E') and not hid.startswith('Ewo'):
				key = 'E'
			else:
				key = hid[:3]
			prefix = metadata.ix[key,'NewPrefix']
			counter[key] += 1
			newid = prefix + str(counter[key]).zfill(3)
			mysites = ' '.join([str(x) for x in sites[i]])
			f.write('%s,%s,%s\n' % (newid, hid, mysites))