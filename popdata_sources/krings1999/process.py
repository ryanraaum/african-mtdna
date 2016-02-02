from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from oldowan.fasta import fasta
from string import translate
import pandas as pd
import sys
import re

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0, 'SeqRange'])

ff = fasta('krings1999.fasta', 'r')
data = ff.readentries()
ff.close()

hids = []
sites = []

for entry in data:
	hids.append(entry['name'])
	sites.append(seq2sites(entry['sequence']))

# some of the sequences are short. Many are just missing a base or two
# from the beginning or the end - will keep those
# a few are missing large chunks of the end, so will drop those
ok = [41, 42, 43, 44, 54, 55, 83, 85, 89, 136, 137, 157, 178]
skip = [8, 123, 124]

# validate
passed_validation = True

for i in range(len(sites)):
	if i in ok or i in skip:
		pass
	else:
		seq1 = data[i]['sequence'].upper()
		if not seq1 == translate(sites2seq(sites[i], region), None, '-'):
				passed_validation = False
				print i, hids[i]

counter = {}
for k in metadata.GroupName:
	counter[k] = 0

if passed_validation:
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			if i not in skip:
				hid = hids[i]
				key = None
				prefix = None
				for pattern in metadata.index:
					if re.match(pattern, hid) is not None:
						key = metadata.ix[pattern, 'GroupName']
						prefix = metadata.ix[pattern, 'NewPrefix']
				counter[key] += 1
				newid = prefix + str(counter[key]).zfill(3)
				mysites = ' '.join([str(x) for x in sites[i]])
				f.write('%s,%s,%s\n' % (newid, hid, mysites))