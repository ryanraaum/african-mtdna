# go through the genbank file, extract ids and population attributions
# drop the mandenka and yoruba (already in db from Lippold et al. 2014)
# and write out in fasta
import sys
from oldowan.fasta import fasta

sys.path.append('../../scripts')
from utils import *
from genbank import read_genbank

entries = read_genbank('barbieri2012.gb', what='filename')

def getpop(x):
	source = x['features'][0][1]
	for x in source:
		if isinstance(x, tuple):
			if x[0] == 'pop_variant':
				return x[1]

def getid(x):
	words = x['definition'].split()
	return words[3]

ff = fasta('barbieri2012.fasta', 'w')

for i in range(len(entries)):
	hid = getid(entries[i])
	if not hid.startswith('MAN') and not hid.startswith('YOR'):
		pop = getpop(entries[i])
		seq = entries[i]['sequence'].upper()
		newentry = {'name': "%s %s" % (hid, pop), 'sequence':seq}
		ff.write(newentry)
		print hid, pop

ff.close()