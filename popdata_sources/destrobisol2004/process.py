from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys
import csv

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

header = []
hids = []
data = []
counts = []
sites = []

with open('destrobisol_tableA1.csv', 'rU') as f:
	reader = csv.reader(f)
	for i in range(4):
		header.append(reader.next())
	for row in reader:
		hids.append(row[0])
		data.append(row[1])
		counts.append(row[2:4])

sitenums = []
for i in range(len(header[0][1])):
	sitenums.append(header[0][1][i] + header[1][1][i] + header[2][1][i])
sitenums = [x for x in sitenums if x != '   ']

for i in range(len(data)):
	values = data[i].split()
	assert len(values) == len(sitenums)

	variants = []
	for j in range(len(values)):
		if values[j] != '.':
			variants.append(sitenums[j] + values[j])

	sites.append(str2sites(' '.join(variants), add16k=True))

newcounts = []
for i in range(len(counts)):
	fixed = []
	for j in range(len(counts[i])):
		if counts[i][j] == '\xc9':
			fixed.append(0)
		else:
			fixed.append(int(counts[i][j]))
	newcounts.append(fixed)

header2 = []
hids2 = []
data2 = []
counts2 = []
sites2 = []

with open('destrobisol_tableA2.csv', 'rU') as f:
	reader = csv.reader(f)
	for i in range(4):
		header2.append(reader.next())
	for row in reader:
		hids2.append(row[0])
		data2.append(row[1])
		counts2.append(row[2])

sitenums2 = []
for i in range(len(header2[0][1])):
	sitenums2.append(header2[0][1][i] + header2[1][1][i] + header2[2][1][i])
sitenums2 = [x for x in sitenums2 if x != '   ']

for i in range(len(data2)):
	values = data2[i].split()
	assert len(values) == len(sitenums2)

	variants = []
	for j in range(len(values)):
		if values[j] != '.':
			variants.append(sitenums2[j] + values[j])

	sites2.append(str2sites(' '.join(variants), add16k=True))

newcounts2 = [int(x) for x in counts2]

## Validate variant sites
passed_validation = True

for i in range(len(sites)):
	seq = sites2seq(sites[i], region)
	mysites = seq2sites(seq)
	if not sites[i] == mysites:
		if not translate(seq, None, '-') == translate(sites2seq(mysites, region), None, '-'):
			passed_validation = False
			print i

for i in range(len(sites2)):
	seq = sites2seq(sites2[i], region)
	mysites = seq2sites(seq)
	if not sites2[i] == mysites:
		if not translate(seq, None, '-') == translate(sites2seq(mysites, region), None, '-'):
			passed_validation = False
			print i

if passed_validation:
	counter = [0,0]
	with open('processed.csv', 'w') as f:
		for i in range(len(hids)):
			seq = sites2seq(sites[i], region)
			seq = translate(seq, None, '-')
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			for j in range(2):
				for k in range(newcounts[i][j]):
					counter[j] += 1
					prefix = metadata.ix[j,'NewPrefix']
					newid = prefix + str(counter[j]).zfill(3)
					f.write('%s,%s,%s\n' % (newid, hids[i], mysites))

		counter = 0
		for i in range(len(hids2)):
			seq = sites2seq(sites2[i], region)
			seq = translate(seq, None, '-')
			mysites = seq2sites(seq)
			mysites = ' '.join([str(x) for x in mysites])
			for j in range(newcounts2[i]):
				counter += 1
				prefix = metadata.ix[2,'NewPrefix']
				newid = prefix + str(counter).zfill(3)
				f.write('%s,%s,%s\n' % (newid, hids2[i], mysites))
