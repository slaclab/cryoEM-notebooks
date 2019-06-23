#import diffmaps
import h5py
import glob
import pylab as pyl
import numpy as np
from scipy.linalg import svd
#from scipy.misc import imresize
from scipy.spatial.distance import pdist, squareform
#
def dm_affinity_matrix(data, b, sigma2, alph):
	"""
	 dm_affinity_matrix
	"""
	n, D = np.shape(data)
	#
	xdist2  = squareform(pdist(data))
	#    
	ids     = np.argsort(xdist2, axis=0)
	xdsts2  = np.sort(xdist2, axis=0)  # sort columns of distance matrix
	xdsts2  = xdsts2[0:b,:]          # keep only `b` nearest neighbors
	ids     = ids[0:b, :]
	# build weight matrix and transition probability matrix
	K = np.zeros((n, n))
	for i in range(n):
		for j in range(b):
			K[i, ids[j,i]] = np.exp( -xdsts2[j,i]**2/sigma2 )

	K = np.maximum(K, K.transpose()) # make it symmetric
	logK = np.log(np.sum(K.flatten())) # for diagnosing epsilon

	Q = np.power(np.sum(K, axis=1), alph)
	W = K/Q/Q[:,np.newaxis]

	Q = np.sum(W, axis=1)
	W = W/Q

	return W
