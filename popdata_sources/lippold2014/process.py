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

## load sample info
sinfo = pd.read_csv('HGDP_info.csv', index_col=0)
newindices = ['HGDP' + str(x).zfill(5) for x in sinfo.index]
sinfo['hgdpid'] = newindices
sinfo = sinfo.set_index('hgdpid')

ff = fasta('hgdp_africa.fasta', 'r')
data = ff.readentries()
ff.close()

hids = []
sites = []

for entry in data:
	words = entry['name'].split()
	hids.append(words[4])
	sites.append(seq2sites(entry['sequence']))

# three sequences have an 'N' at around 309 that breaks validation
# this will be treated as a heteroplasy of the T there and ignored
skip = [64, 67, 73]

# validate
passed_validation = True

for i in range(len(sites)):
	if i not in skip:
		seq1 = data[i]['sequence'].upper()
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
			key = sinfo.ix[hid,'PopulationName']
			prefix = metadata.ix[key,'NewPrefix']
			counter[key] += 1
			newid = prefix + str(counter[key]).zfill(3)
			mysites = ' '.join([str(x) for x in sites[i]])
			f.write('%s,%s,%s\n' % (newid, hid, mysites))