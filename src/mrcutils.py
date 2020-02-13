import shutil, os, sys
#from prody import *
from pylab import *
from matplotlib import pyplot as plt
import numpy as np
import mrcfile
from scipy.ndimage import gaussian_filter

# routines on manipulation of mrc data

def mrc_stats(mrc_filename, get='std'):
    """ mrcstats
    """
    data = mrc2data(mrc_filename)
    mean = np.mean(data)
    if(get=='mean'):
        value = mean
    elif(get=='std'):
        value = np.std(data-mean)
    elif(get=='min'):
        value = np.min(data)
    elif(get=='max'):
        value = np.max(data)
    elif(get=='sum'):
        value = np.sum(data)
    return value

def mrc_algebra(mrc1,mrc2,mrc_out,operation='add'):
    """mrc_algebra: mrc_out = mrc1 operation mrc2
    """
    data1 = mrc2data(mrc1)
    data2 = mrc2data(mrc2)
    if(operation=='add'):
        data = data1 + data2
    elif(operation=='subtract'):
        data = data1 - data2
    data2mrc(mrc_out,data,mrc_template=mrc1)
        
def mrc_select(mrc_filename, mode='above_value', value=0.):
    """mrc_select
    """
    data = mrc2data(mrc_filename)
    if(mode=='above_value'):
        data_selected = np.where(data >  value, data, 0.0)
    elif(mode=='equal_value'):
        data_selected = np.where(data == value, data, 0.0)
    else:
        data_selected = np.where(data <  value, data, 0.0)
    return data_selected
        
# routines mrc2bla or bla2mrc
    
def mrc2mask(mrc_filename, mask_filename, sigma_blur=0., threshold=0.1):
    """ mrc2mask: set to 1 any non-zero value, blurs, and binarize around threshold
    """
    data = mrc2data(mrc_filename)
    mask = data2mask(data, sigma_blur=sigma_blur, threshold=threshold)
    data2mrc(mask_filename, mask, mrc_template=mrc_filename)

def seg2mask(input_seg, output_key, sigma_blur = 0., threshold=0.1,sort=None,verbose=False):
    """seg2mask
    """
    segments = mrc2data(input_seg)
    domains  = np.unique(segments).astype(int)
    index_order = domains
    if sort is not None:
        if(sort=='volume'):
            volume = []
            if (domains.shape[0] > 1):
                for i in domains:
                    data_domain = mrc_select(input_seg, mode='equal_value', value=i)
                    masked = np.ma.masked_equal(data_domain, 0)
                    volume.append(np.sum(masked))
            index_order = np.argsort(volume)[::-1]
    if (domains.shape[0] > 1):
        for i in domains:
            if(i>0):
                data_domain = mrc_select(input_seg, mode='equal_value', value=index_order[i])
                mask = data2mask(data_domain, sigma_blur=sigma_blur, threshold=threshold)
                data2mrc(output_key+str(i)+'.mrc', mask, mrc_template=input_seg)
                if verbose:
                    print("{0} > volume = {1}".format(output_key+str(i)+'.mrc',volume[index_order[i]]))
                
def seg2kept(input_seg, output_key, keep_mask):
    """seg2kep
    """
    ikeep=0
    segments = mrc2data(input_seg)
    domains  = np.unique(segments).astype(int)
    if (domains.shape[0] > 1):
        for i in domains:
            if(i>0):
                if keep_mask[i-1]:
                    ikeep += 1
                    shutil.copy(output_key+str(i)+'.mrc',output_key+'kept_'+str(ikeep)+'.mrc')
                os.remove(output_key+str(i)+'.mrc')
    print("keeping ",ikeep," domains")
    return ikeep
                
def data2mask(data, sigma_blur=0., threshold=0.1):
    """data2mask: set to 1 any non-zero value, blurs, and binarize around threshold
    """
    mask = np.where(data > 0, 1.0, 0.0)
    if(sigma_blur > 0):
        mask = gaussian_filter(mask, sigma_blur)
        mask = np.where(mask > threshold, 1.0 , 0.0 ) #.astype(np.int8)
    return mask.astype(np.int8)
    
def mrc2data(mrc_filename):
    """ mrc2data
    """
    mrc  = mrcfile.open(mrc_filename, mode='r+')
    data = mrc.data
    mrc.close()
    return data

def data2mrc(mrc_filename,data,mrc_template=None):
    """ data2mrc
    """
    if mrc_template is None:
        print("Please provide a .mrc template to update data from")
    else:
        shutil.copy(mrc_template, mrc_filename)
        mrc = mrcfile.open(mrc_filename, mode='r+')
        mrc.data[:] = data
        mrc.close()

### mrc and their power spectra

def data2psd(data, log=True):
    """data2psd
    """
    data_psd = np.fft.fftn(data, s=data.shape)
    data_psd = np.fft.fftshift(data_psd)
    data_psd = np.absolute(data_psd)**2
    if log:
        data_psd = np.log(data_psd)
    return data_psd

def showdata(data, level=[1,1]):
    """
    """
    vmin = np.mean(data) - level[0]*np.std(data)
    vmax = np.mean(data) + level[1]*np.std(data)
    fig = plt.figure(figsize=(12,4))
    plt.subplot(1,3,1)
    plt.imshow(np.mean(data, axis=0), vmin=vmin, vmax=vmax, cmap='Greys_r')
    plt.subplot(1,3,2)
    plt.imshow(np.mean(data, axis=1), vmin=vmin, vmax=vmax, cmap='Greys_r')
    plt.subplot(1,3,3)
    plt.imshow(np.mean(data, axis=2), vmin=vmin, vmax=vmax, cmap='Greys_r')
    plt.show()

def data2radialprofile(data):
    """
    adapted from https://gist.github.com/ViggieSmalls/3bc5ec52774cf6e672f49723f0aa4a47
    """
    center = (np.floor(data.shape)/2).astype(int)
    indices = np.indices(data.shape)
    indices_centered = (indices.T - center).T
    radius = np.linalg.norm(indices_centered, axis=0)
    radius_indices = np.argsort(radius.flat)
    radius_sorted  = radius.flat[radius_indices]
    data_sorted    = data.flat[radius_indices]
    radius_integer = radius_sorted.astype(np.int16)
    deltar = radius_integer[1:] - radius_integer[:-1]  # assume all radii represented (more work if not)
    rind = np.where(deltar)[0]  # location of changed radius
    nr = rind[1:] - rind[:-1]  # number of pixels in radius bin
    data_cumul = np.cumsum(data_sorted, dtype=np.float64)  # cumulative sum for increasing radius
    tbin = data_cumul[rind[1:]] - data_cumul[rind[:-1]]  # sum for image values in radius bins
    radialprofile = tbin / nr
    return radialprofile


### ad-hoc manipulations

def cut_thresholded_map(input_map, 
                        low_threshold=1, low_blur=1, low_bin=0.1, 
                        high_threshold=6, high_blur=15, high_bin=0.01):
    """
    """
    data_dry = mrc_select(input_map, mode='above_value', value=high_threshold)
    mask_dry = data2mask(data_dry, sigma_blur=high_blur, threshold=high_bin)
    data_fat = mrc_select(input_map, mode='above_value', value=low_threshold)
    mask_fat = data2mask(data_fat, sigma_blur=low_blur, threshold=low_bin)
    body0 = np.minimum(mask_dry,mask_fat)
    bodyK = mask_fat - body0
    return body0, bodyK

### visualize

def view_map_xyzproj(data,style=None):
    masked = np.ma.masked_equal(data, 0)
    N = data.shape[0]*data.shape[1]*data.shape[2]
    K = int(data.shape[0]/5)
    if style is None:
        fig = plt.figure(figsize=(12,8))
    else:
        fig = plt.figure(figsize=(12,4))
    if style is None:
        plt.subplot(231)
    else:
        plt.subplot(131)
    plt.imshow(np.mean(masked,axis=0))
    plt.colorbar()
    if style is None:
        plt.subplot(232)
    else:
        plt.subplot(132)
    plt.imshow(np.mean(masked,axis=1))
    plt.colorbar()
    if style is None:
        plt.subplot(233)
    else:
        plt.subplot(133)
    plt.imshow(np.mean(masked,axis=2))
    plt.colorbar()
    if style is None:
        plt.subplot(212)
        plt.hist(masked.reshape(N),bins=K,log=True,histtype='step')
    plt.tight_layout()
