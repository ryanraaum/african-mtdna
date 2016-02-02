from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from oldowan.polymorphism import Polymorphism
from string import translate
import pandas as pd
import numpy as np
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

region_str = metadata.ix[0,'SeqRange']
region_parts = region_str.split(';')
region1 = range2region(region_parts[0])
region2 = range2region(region_parts[1])

with open('chen2000.txt', 'rU') as f:
	data = f.readlines()

iids = []
counts = np.zeros((len(data), 2), dtype=np.int)
sites = []

for l in data:
	parts = l.strip().split()
	iids.append(parts[0])
	row = len(iids) - 1
	counts[row,] = [int(x) for x in parts[1:3]]
	sites.append(' '.join(parts[3:]))

## Validate
passed_validation = True

# there are sites in the source table that are not actual variant sites
# sequence 9 (index 8) has 263A as a variant
# sequence 12 (index 11) has 16223C as a variant
not_polys = [Polymorphism(263,0,'A'), Polymorphism(16223,0,'C')]

for i in range(len(sites)):
	curr_sites = sites[i]
	curr_polys = [x for x in str2sites(curr_sites) if x not in not_polys]
	cseq1 = sites2seq(curr_sites, region1)
	cseq2 = sites2seq(curr_sites, region2)
	mysites1 = seq2sites(cseq1)
	mysites2 = seq2sites(cseq2)
	mysites = mysites1 + mysites2
	if not mysites == curr_polys:
		passed_validation = False
		print iids[i]

if passed_validation:
	counters = [1] * 2
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			curr_sites = str2sites(sites[i])
			mysites = [x for x in curr_sites if x not in not_polys]
			mysites = ' '.join([str(x) for x in mysites])
			for j in range(2):
				if counts[i,j] > 0:
					prefix = metadata.ix[j,'NewPrefix']
					for k in range(counts[i,j]):
						newid = prefix + str(counters[j]).zfill(3)
						f.write('%s,%s,%s\n' % (newid, iids[i], mysites))
						counters[j] += 1