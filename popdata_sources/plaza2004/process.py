from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region('16024-16569;1-322')

locations = []
hids = []
sites = []

with open('plaza2004.txt', 'r') as f:
	line = f.readline()
	parts = line.split()
	hvr1_1 = parts[0]
	hvr2_1 = parts[1]

	line = f.readline()
	parts = line.split()
	hvr1_2 = parts[0]
	hvr2_2 = parts[1]

	line = f.readline()
	parts = line.split()
	hvr1_3 = parts[0]
	hvr2_3 = parts[1]

	line = f.readline()
	parts = line.split()
	hvr1_4 = parts[0]

	line = f.readline()
	parts = line.split()
	hvr1_5 = parts[0]

	for i in range(len(hvr1_1)):
		locations.append("%s%s%s%s%s" % (hvr1_1[i], hvr1_2[i], hvr1_3[i], hvr1_4[i], hvr1_5[i]))

	for i in range(len(hvr2_1)):
		locations.append("%s%s%s" % (hvr2_1[i], hvr2_2[i], hvr2_3[i]))

	locations[-2] = locations[-2] + '.1'
	locations[-3] = locations[-3] + '.2'
	locations[-4] = locations[-4] + '.1'

	f.readline() # skip past anderson sequence

	for line in f:
		parts = line.split()
		# drop non-Mbundu (AN9) and individuals without HVR2 (AN130, AN42)
		if parts[0] not in ['AN9', 'AN130', 'AN42']:
			hids.append(parts[0])
			bits = [x for x in parts[1]] + [x for x in parts[2]]
			assert len(bits) == len(locations)
			variants = []
			for i in range(len(locations)):
				if bits[i] != '.':
					variants.append(locations[i]+bits[i])
			sites.append(str2sites(' '.join(variants)))

## Validate
passed_validation = True

for i in range(len(hids)):
	seq = sites2seq(sites[i], region)
	mysites = seq2sites(seq)
	if not sites[i] == mysites:
		if not translate(seq, None, '-') == translate(sites2seq(mysites, region), None, '-'):
			passed_validation = False
			print i, hids[i]


if passed_validation:
	counter = 0
	prefix = metadata.ix[0,'NewPrefix']
	with open('processed.csv', 'w') as f:
		for i in range(len(hids)):
			counter = counter + 1
			newid = prefix + str(counter).zfill(3)
			seq = sites2seq(sites[i], region)
			seq = translate(seq, None, '-')
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			f.write('%s,%s,%s\n' % (newid, hids[i], mysites))