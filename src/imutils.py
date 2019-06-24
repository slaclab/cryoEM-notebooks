# 
# all utilities to handle images: binning, normalization, ...
#
import numpy as np
import skimage
#
def binstack(imstack, bin_size=2, method='skimage'):
	""" binning: provides several approaches to binning images in a stack
	
	INPUTS:
	======
	- imstack  : 
	first axis gives image index in stack.
	- bin_size : 
	reducing factor in one direction. 
	- method   : 
	. 'skimage' 
	"""
	L, Nrow, Ncol = imstack.shape
	# make sure the new image size is a multiple of 2
	newNrow = 2*int(np.ceil((Nrow/bin_size)/2))
	newNcol = 2*int(np.ceil((Ncol/bin_size)/2))
	#
	if(method=='skimage'):
		binned_imstack = np.zeros((imstack.shape[0],newNrow,newNcol))
		for i in np.arange(imstack.shape[0]):
			binned_imstack[i,...] = skimage.transform.resize(imstack[i,...], 
                                                                         output_shape=(newNrow,newNcol), 
                                                                         order=1,
                                                                         mode='edge',
                                                                         preserve_range=True,
                                                                         anti_aliasing=True)
	return binned_imstack
#
