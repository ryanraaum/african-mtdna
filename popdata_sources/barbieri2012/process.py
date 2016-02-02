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

ff = fasta('barbieri2012.fasta', 'r')
data = ff.readentries()
ff.close()

hids = []
groups = []
sites = []
sequences = []

for entry in data:
	words = entry['name'].split()
	if entry['sequence'].count('N') > 10:
		print "Too many Ns in %s (%s), skipping" % (words[0], words[1])
	else:
		hids.append(words[0])
		groups.append(words[1])
		sites.append(seq2sites(entry['sequence']))
		sequences.append(entry['sequence'])

# validate
passed_validation = True

for i in range(len(sites)):
	seq1 = sequences[i]
	seq2 = translate(sites2seq(sites[i], region), None, '-')
	if not seq1 == seq2:	
		# some sequences have an N in the poly-C region around 309.
		if seq1[:305] == seq2[:305] and seq1[-16254:] == seq2[-16254:]:
			pass
		else:
			passed_validation = False
			print i, hids[i]

counter = {}
for k in metadata.index:
	counter[k] = 0

if passed_validation:
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			key = groups[i]
			prefix = metadata.ix[key,'NewPrefix']
			counter[key] += 1
			newid = prefix + str(counter[key]).zfill(3)
			mysites = ' '.join([str(x) for x in sites[i]])
			f.write('%s,%s,%s\n' % (newid, hid, mysites))