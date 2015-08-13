from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from oldowan.fasta import fasta
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
region = range2region(metadata.ix[0,'SeqRange'])

groups = pd.read_csv('barbieri2014_groups.csv', index_col=0)
newindex = [x.upper() for x in groups.index]
groups.index = pd.Index(newindex)

hids = []
seqs = []

for filename in ['barbieri2013a.fasta', 'barbieri2013b.fasta', 'barbieri2013c.fasta', 'barbieri2014.fasta']:
	ff = fasta(filename, 'r')
	data = ff.readentries()
	ff.close()

	for e in data:
		if e['sequence'].count('N') < 20:
			hid = e['name'].split()[4].upper()
			if hid in groups.index:
				hids.append(hid)
				seqs.append(e['sequence'])

sites = []
for s in seqs:
	sites.append(seq2sites(s, ambig_cutoff=20))

## Validate
passed_validation = True

for i in range(len(sites)):
	hid = hids[i]
	seq = translate(sites2seq(sites[i], region), None, '-')
	if not seq == seqs[i]:
		# some sequences have N in position 308 and this doesn't parse properly in my converter
		# not going to worry about it because will always drop the variants in this region
		if not seq[:305] == seqs[i][:305] and seq[-16260:] == seqs[i][-16260:]:
			passed_validation = False
			print i, hids[i]

if passed_validation:
	counter = {}
	for k in metadata.index:
		counter[k] = 0
	with open('processed.csv', 'w') as f:
		for i in range(len(sites)):
			hid = hids[i]
			grp = groups.ix[hid,'Assigned population']
			mysites = ' '.join([str(x) for x in sites[i]])
			prefix = metadata.ix[grp,'NewPrefix']
			counter[grp] += 1
			num = str(counter[grp]).zfill(3)
			newid = prefix + num
			f.write('%s,%s,%s\n' % (newid, hid, mysites))