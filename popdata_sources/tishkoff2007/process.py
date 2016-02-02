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

ff = fasta('tishkoff2007.fasta', 'rU')
data = ff.readentries()
ff.close()

drop = ['San_43', 'San_67', 'tzbg040', 'tzdt045', 'tzhz108', 'tzhz130', 'tzhz131']

hids = []
seqs = []

for e in data:
	name_parts = e['name'].split()
	if name_parts[0] not in drop:
		hids.append(name_parts[0])
		seqs.append(e['sequence'])

## Validate
passed_validation = True

for i in range(len(seqs)):
	mysites = seq2sites(seqs[i])
	myseq = translate(sites2seq(mysites, region), None, '-')
	if not seqs[i] == myseq:
		passed_validation = False
		print i, hids[i]

if passed_validation:
	with open('processed.csv', 'w') as f:
		for i in range(len(seqs)):
			mysites = ' '.join([str(x) for x in seq2sites(seqs[i])])
			origid = hids[i]
			prefix = metadata.ix[origid[:4],'NewPrefix']
			num = origid[4:].split('_')[0].zfill(3)
			newid = prefix + num
			f.write('%s,%s,%s\n' % (newid, origid, mysites))