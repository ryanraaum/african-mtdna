from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys
import csv

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region('16051-16569;1-340')

hvr1 = []
hvr2 = []
counts = []
sites = []

with open('trovoada2004.csv', 'rU') as f:
	reader = csv.reader(f)
	reader.next() # skip past header
	for row in reader:
		counts.append(row[:3])
		hvr1.append(row[3].split())
		hvr2.append(row[4].split())

# fix counts
def fix(x):
	if x == ' ':
		return 0
	return int(x)

newcounts = []
for i in range(len(counts)):	
	newcounts.append([fix(x) for x in counts[i]])

# transversions are reported as A/D (rcrs/derived), which is not standard notation
# also, C insertions at 309 and 315 are reported as 309.1 and 315.1, not 309.1C and 315.1C
# so need to go through and fix all of those
# also, sites outside of the stated sequence range are in parentheses, i.e. (385) - drop those
def fix2(x):
	if '/' in x:
		derived = x[-1]
		x = x[:-3] + x[-1]
	elif x == '309.1':
		x = '309.1C'
	elif x == '315.1':
		x = '315.1C'
	elif '(' in x:
		x = ''
	return x	

for i in range(len(hvr1)):
	hvr1[i] = [fix2(x) for x in hvr1[i]]

for i in range(len(hvr2)):
	hvr2[i] = [fix2(x) for x in hvr2[i]]

for i in range(len(hvr1)):
	sites.append(str2sites(' '.join(hvr1[i]), add16k=True) + str2sites(' '.join(hvr2[i])))

## Validate variant sites
passed_validation = True

for i in range(len(sites)):
	seq = sites2seq(sites[i], region)
	mysites = seq2sites(seq)
	if not sites[i] == mysites:
		if not translate(seq, None, '-') == translate(sites2seq(mysites, region), None, '-'):
			passed_validation = False
			print i

if passed_validation:
	counter = [0] * 3
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = "H%d" % (i + 1,)
			seq = sites2seq(sites[i], region)
			seq = translate(seq, None, '-')
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			for j in range(3):
				for k in range(newcounts[i][j]):
					counter[j] += 1
					newid = metadata.ix[j,'NewPrefix'] + str(counter[j]).zfill(3)
					f.write('%s,%s,%s\n' % (newid, hid, mysites))