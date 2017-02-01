import numpy as np
from scipy.io import loadmat  # this is the SciPy module that loads mat-files
import pandas as pd

def load(infile, dataset=None, idKey='id'):
	mat = loadmat(infile)  # load mat-file
	#print mat.keys()
	if dataset == None:
		for key in mat.keys():
			if key[0] != '_':
				dataset = key
				break
		print 'Assuming you mean to load',dataset,'from',infile
	mdata = mat[dataset]  # variable in mat file
	mdtype = mdata.dtype  # dtypes of structures are "unsized objects"
	ndata = {n: mdata[n][0, 0] for n in mdtype.names}

	# Reconstruct the columns of the data table from just the time series
	# Use the number of intervals to test if a field is a column or metadata
	columns = [n for n, v in ndata.iteritems()]

	# There were supposed to be useful functions instead of this but those 'complained'
	# going the rude but working way now
	dataDict = dict()
	for c in columns:
		curCol = []

		for thing in mdata:
			value = thing[c][0][0]
			if type(value) == np.unicode_:
				value = str(value)
			if c == 'vp_chr':
				value = value[0]
			curCol.append(value)
		dataDict[c] = curCol

	df = pd.DataFrame(dataDict,
		index=dataDict[idKey],
		columns=columns)

	return df

def listIndexes(args):
	df = load(args.infile)
	for item in df.index:
		print str(item)

#df = load('/mnt/hpc/Dataset_info.mat')
#print df.loc['Upp2']
