import numpy as np
import numpy.ma as ma
import cPickle as pickle
import gzip

# Based on code posted at http://werthmuller.org/blog/2014/basemap/
# Load the population density file to get header information
popdensname = 'misc/afds00ag.asc'
popdens_file = open(popdensname, 'r')

# Read header (number of columns and rows, cell-size, and lower left coordinates)
ncols = int(popdens_file.readline().split()[1])
nrows = int(popdens_file.readline().split()[1])
xllcorner = float(popdens_file.readline().split()[1])
yllcorner = float(popdens_file.readline().split()[1])
cellsize = float(popdens_file.readline().split()[1])
nodata = float(popdens_file.readline().split()[1])
popdens_file.close()

# Read in population density data as a whole, disregarding first five rows (header)
popdens = np.loadtxt(popdensname, skiprows=6)

# Swap the rows
popdens[:nrows+1, :] = popdens[nrows+1::-1, :]
popdens[popdens == -9999] = np.nan

# Create longitude and latitude vectors for popdens
lons = np.arange(xllcorner, xllcorner+cellsize*ncols, cellsize)
lats = np.arange(yllcorner, yllcorner+cellsize*nrows, cellsize)

# mask out invalid data
mpopdens = ma.masked_invalid(popdens)

# combine it all togther into a dict
pdd = {'lons':lons, 'lats':lats, 'popdens':mpopdens}

with gzip.open('misc/popdens.pkl.gz', 'wb') as f:
	pickle.dump(pdd, f, protocol=2)
