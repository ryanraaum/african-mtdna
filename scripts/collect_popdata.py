import os
from shutil import copyfile
import pandas as pd
import numpy as np

##### Collect data from subdirs

subdirs = [name for name in os.listdir('popdata_sources') 
			if os.path.isdir(os.path.join('popdata_sources', name))]

df = open('popdata/data.csv', 'w')
mdf = open('popdata/metadata.csv', 'w')

firstmdf = True

for d in subdirs:
	d = os.path.join('popdata_sources', d)
	with open(os.path.join(d, 'processed.csv'), 'rU') as f:
		for l in f:
			df.write(l)

	with open(os.path.join(d, 'metadata.csv'), 'rU') as f:
		if firstmdf:
			firstmdf = False
		else:
			f.readline() # skip past header
		for l in f:
			mdf.write(l)
		mdf.write('\n')

df.close()
mdf.close()

### Update metadata

data = pd.read_csv('popdata/data.csv', index_col=0)
prefixes = pd.Series([x[:-3] for x in data.index])
counts = prefixes.value_counts()

metadata = pd.read_csv('popdata/metadata.csv', index_col=1)
metadata = metadata.drop('Selector', 1)
#metadata['N'] = pd.Series(np.zeros(metadata.shape[0]), index=metadata.index, dtype='int64')
N = pd.Series(np.zeros(metadata.shape[0]), index=metadata.index, dtype='int64')

for i in counts.index:
	#metadata.ix[i,'N'] = counts[i]
	N[i] = counts[i]

metadata.insert(1, 'N', N)

metadata.to_csv('popdata/metadata.csv')

