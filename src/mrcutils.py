"""
MRCUTILS - F.Poitevin, Stanford, 2019
a collection of python tools to manipulate images
"""
import matplotlib
from matplotlib import pyplot as plt
import math
import numpy as np
import scipy.misc
from scipy import ndimage
import cairo
from skimage import img_as_float
import mrcfile
import hyperspy.api as hs
hs.preferences.GUIs.warn_if_guis_are_missing = False
hs.preferences.save()
#
def mrc2array(mrcfile,squarify=False,n_stack=256,n_ovl=0):
    """ mrc2array: reads a mrc file and returns image as numpy array
    """
    data_mrc = hs.load(mrcfile, lazy=True)
    data_array = np.asarray(data_mrc.data[0,...])
    print("Image size: ",data_array.shape)
    if squarify:
        ny, nx = data_array.shape
        n_split = math.floor((np.minimum(nx,ny)-2*n_ovl)/(n_stack))
        n = 2*n_ovl + n_split*n_stack
        data_array = data_array[0:n,0:n]
    print("Image size: ",data_array.shape)
    return data_array

#
def slicenstack(data, n_stack=256, n_ovl=0):
    """ slicenstack: reads a 2D numpy array and returns a 3D numpy array
    """
    if n_ovl == 0:
        data_stack = blockshaped(data, n_stack, n_stack)
    else:
        n_split = math.floor((data.shape[0]-2*n_ovl)/(n_stack))
        n_dilat = n_stack + 2*n_ovl
        data_stack = np.zeros((n_split*n_split,n_dilat,n_dilat))
        print("Array dimensions: ",data_stack.shape)
        i_stack = 0
        for i in np.arange(n_split):
            for j in np.arange(n_split):
                istart = i*n_stack
                istop  = istart + n_dilat
                jstart = j*n_stack
                jstop  = jstart + n_dilat
                rows    = np.arange(istart,istop)
                columns = np.arange(jstart,jstop)
                data_tmp = data[np.ix_(rows,columns)]
                data_stack[i_stack,...] = data_tmp[np.newaxis,...] 
                i_stack += 1
    print("Array dimensions: ",data_stack.shape)
    return data_stack

def unstack(stack,n_ovl=0):
    """ unstack: reads a 3D numpy array and returns a 2D numpy array
    """
    n_split = int(np.sqrt(stack.shape[0]))
    if n_ovl == 0:
        n_stack = stack.shape[1]
        n = n_split*n_stack
        data_unstacked = unblockshaped(stack, n, n)
    else:
        n_stack = stack.shape[1] - 2*n_ovl
        n = n_split*n_stack
        stacks  = np.arange(0,stack.shape[0])
        rows    = np.arange(n_ovl,n_ovl+n_stack)
        columns = np.arange(n_ovl,n_ovl+n_stack)
        stack_new = stack[np.ix_(stacks,rows,columns)]
        print(stack_new.shape)
        data_unstacked = unblockshaped(stack_new, n, n)
    print("Array dimensions: ",data_unstacked.shape)
    return data_unstacked

def update_data(mrcfile,npy_data):
    """ update_data: replace data in mrcfile with npy_data. Returns new mrc image
    """
    image = hs.load(mrcfile, lazy=True)
    image.data = npy_data[np.newaxis,...]
    return image

def psd(data):
    """
    """
    data_psd = np.fft.fft2(data)
    data_psd = np.fft.fftshift(data_psd)
    data_psd = np.log(np.absolute(data_psd)**2)
    return data_psd

#
def prepare_ctffind(input_file,mrcin='micrograph.mrc',mrcout='ctf.mrc',psize=1.0,voltage=200.0,Cs=2.70,Ac=0.07,size=512,minres=50.0,maxres=4.0,mindef=5000,maxdef=50000,stepdef=100):
    """
    """
    f = open(input_file,"w")
    f.write("%s\n" % mrcin)
    f.write("%s\n" % mrcout)
    f.write("%s\n" % str(psize))
    f.write("%s\n" % str(voltage))
    f.write("%s\n" % str(Cs))
    f.write("%s\n" % str(Ac))
    f.write("%s\n" % str(size))
    f.write("%s\n" % str(minres))
    f.write("%s\n" % str(maxres))
    f.write("%s\n" % str(mindef))
    f.write("%s\n" % str(maxdef))
    f.write("%s\n" % str(stepdef))
    f.write("no\nno\nno\nno\nno\n")
    f.close()
# 
def lowpass(data,n):
    """ lowpass: reads a 2D numpy array and returns array of same size
    """
    kernel = np.ones((n,n))
    return ndimage.convolve(data, kernel)

# Found how to slice image in a stack of smaller images here:
# https://stackoverflow.com/questions/16856788/slice-2d-array-into-smaller-2d-arrays
# https://stackoverflow.com/questions/42297115/numpy-split-cube-into-cubes/42298440#42298440

def blockshaped(arr, nrows, ncols):
    """
    Return an array of shape (n, nrows, ncols) where
    n * nrows * ncols = arr.size

    If arr is a 2D array, the returned array looks like n subblocks with
    each subblock preserving the "physical" layout of arr.
    """
    h, w = arr.shape
    return (arr.reshape(h//nrows, nrows, -1, ncols)
               .swapaxes(1,2)
               .reshape(-1, nrows, ncols))


def unblockshaped(arr, h, w):
    """
    Return an array of shape (h, w) where
    h * w = arr.size

    If arr is of shape (n, nrows, ncols), n sublocks of shape (nrows, ncols),
    then the returned array preserves the "physical" layout of the sublocks.
    """
    n, nrows, ncols = arr.shape
    return (arr.reshape(h//nrows, -1, nrows, ncols)
               .swapaxes(1,2)
               .reshape(h, w))

###################
# SYNTHETIC STUFF #
###################

def create_mrc(mrcfile_out='out.mrc',signal_data=None,noise_data=None,signal_strength=1.0,noise_strength=0.0):
    """ create_mrc:
    """
    if signal_data is None or noise_data is None:
        print("error: provide signal and noise arrays")
    else:
        if signal_data.shape != noise_data.shape:
            print("error: provide signal and noise arrays of same shape")
        else:
            measure = signal_strength*signal_data
            measure += noise_strength*noise_data
            print("... about to write ",mrcfile_out)
            with mrcfile.new(mrcfile_out,overwrite=True) as mrc:
                mrc.set_data(np.zeros(measure.shape, dtype=np.float32))
                mrc.data[:,:] = measure

def draw_in_npy(data):
    # see https://stackoverflow.com/a/10031877/9023182
    nx = data.shape[0]
    ny = data.shape[1]
    n = np.minimum(nx,ny)
    data = np.zeros((data.shape), dtype=np.uint8)
    surface = cairo.ImageSurface.create_for_data(data, cairo.FORMAT_ARGB32, nx, ny)
    cr = cairo.Context(surface)

    # fill with solid black
    cr.set_source_rgb(0.0, 0.0, 0.0)
    cr.paint()

    # draw white circle
    cr.arc(math.floor(n/4), math.floor(n/4), 400, 0, 2*math.pi)
    cr.set_line_width(50)
    cr.set_source_rgb(1.0, 1.0, 1.0)
    cr.stroke()
    
    # draw white circle
    cr.arc(math.floor(3*n/4), math.floor(n/2), 400, 0, 2*math.pi)
    cr.set_line_width(30)
    cr.set_source_rgb(1.0, 1.0, 1.0)
    cr.stroke()
    
    # draw white rectangle
    cr.rectangle(math.floor(n/2), math.floor(n/2),400,800)
    cr.set_line_width(30)
    cr.set_source_rgb(1.0, 1.0, 1.0)
    cr.stroke()
    
    # draw white rectangle
    cr.rectangle(math.floor(n/2)-400, math.floor(3*n/4),800,200)
    cr.set_line_width(30)
    cr.set_source_rgb(1.0, 1.0, 1.0)
    cr.stroke()

    # write output
    #print(data[38:48, 38:48, 0])
    #surface.write_to_png("circle.png")
    return img_as_float(data[:,:,0])

