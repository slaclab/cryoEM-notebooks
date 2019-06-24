#
import numpy as np
from matplotlib import pyplot as plt
#
def view_image(img):
	"""
	"""
	fig = plt.figure(figsize=(6,6))
	plt.imshow(img, cmap='Greys')
	plt.colorbar()
#
def view_particles(data, slicing=(1,1,1), figsize=1, ncol=5 ):
	"""
	"""
	view = data[::slicing[0],::slicing[1],::slicing[2]]
	figsize = int(figsize*ncol)
	nrow = np.ceil(view.shape[0]/ncol)
	fig = plt.figure( figsize=(ncol*figsize,nrow*figsize) )
	for i in np.arange( view.shape[0] ):
		fig.add_subplot( nrow, ncol, i+1 )
		plt.imshow(view[i], cmap='Greys')
	plt.tight_layout()
	plt.show()
