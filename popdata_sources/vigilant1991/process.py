from oldowan.mtconvert import seq2sites, sites2seq, str2sites
from string import translate
import pandas as pd
import sys

sys.path.append('../../scripts')
from utils import *
from genbank import read_genbank

## load metadata
metadata = pd.read_csv('metadata.csv', index_col=0)
regionstr = metadata.ix[0,'SeqRange']
regionparts = regionstr.split(';')
region1 = range2region(regionparts[0])
region2 = range2region(regionparts[1])

counts = {}
with open('vigilant1991_counts.txt', 'rU') as f:
	for line in f:
		parts = line.strip().split(',')
		counts[parts[0]] = int(parts[1])

entries = read_genbank('vigilant1991.txt', what='filename')

hids = []
pops = []
seqs = []
sites = []

for i in range(len(entries)):
	e = entries[i]
	name = e['definition'].split()
	if e['sequence'].count('n') > 10:
		print 'skipping isolate %s' % name[3]
	else:
		hids.append(name[3])
		seqs.append(e['sequence'])
		if name[3] == '63':
			pops.append('Yoruban')
		else:
			for f in e['features']:
				if f[0] == 'source':
					for f2 in f[1]:
						if isinstance(f2, tuple):
							if f2[0] == 'note':
								pops.append(f2[1])
		if 'complete' in e['definition'] or name[3] in ['3', '65']:
			sites.append(seq2sites(e['sequence']))
		else:
			found = False
			for f in e['features']:
				if f[0] == 'misc_feature' and f[1][1][1] == 'segment 1':
					c = int(f[1][0].split('..')[1]) - 1
					s = e['sequence']
					sites.append(seq2sites(s[:c]) + seq2sites(s[c:]))
					found = True
			if not found:
				print 'problem with isolate %s' % name[3]

# Vigilant GenBank data have variable sequence lengths
# normalize all sites to specified range in the metadata file
mysites = []

for i in range(len(sites)):
	seq1 = translate(sites2seq(sites[i], region1), None, '-')
	seq2 = translate(sites2seq(sites[i], region2), None, '-')
	s = seq2sites(seq1) + seq2sites(seq2)
	mysites.append(' '.join([str(x) for x in s]))

counter = {}
for k in metadata.index:
	counter[k] = 0
with open('processed.csv', 'w') as f:
	for i in range(len(sites)):
		hid = hids[i]
		grp = pops[i]
		if grp in metadata.index:
			repeat = 1
			if hid in counts.keys():
				repeat = counts[hid]
			for j in range(repeat):
				prefix = metadata.ix[grp,'NewPrefix']
				counter[grp] += 1
				num = str(counter[grp]).zfill(3)
				newid = prefix + num
				f.write('%s,%s,%s\n' % (newid, hid, mysites[i]))