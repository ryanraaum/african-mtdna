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
region = range2region(metadata.ix['Nai', 'SeqRange'])

ff = fasta('brandstaetter2004.fasta', 'r')
data = ff.readentries()
ff.close()

## Validate
passed_validation = True

for entry in data:
	seq1 = entry['sequence']
	# Brandstatter et al. put an N at the end of an unstable poly-C run at the end
	#    of HVR3 in 5 samples. This spurious N messes with my conversion utility,
	#    so I strip it out.
	if seq1.endswith('NACA'):
		seq1 = seq1[:-4] + 'ACA'
	mysites = seq2sites(seq1)
	if not seq1 == translate(sites2seq(mysites, region), None, '-'):
			passed_validation = False
			print entry['name']

if passed_validation:
	with open('processed.csv', 'w') as f:
		for entry in data:
			name_parts = entry['name'].split()
			origid = name_parts[0]
			key = name_parts[0][:3]
			m = re.search(r'[0-9]+', name_parts[0])
			counter = m.group()
			newid = metadata.ix[key,'NewPrefix'] + counter
			seq = entry['sequence']
			# Brandstatter et al. put an N at the end of an unstable poly-C run at the end
			#    of HVR3 in 5 samples. This spurious N messes with my conversion utility,
			#    so I strip it out.
			if seq.endswith('NACA'):
				seq = seq[:-4] + 'ACA'
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			f.write('%s,%s,%s\n' % (newid, origid, mysites))