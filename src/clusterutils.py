import numpy as np
from matplotlib import pyplot as plt
from sklearn import mixture
import skimage.morphology as morphology
from skimage.morphology import watershed
from skimage.feature import peak_local_max
from scipy import ndimage
from scipy import stats

## BINNING

def bin_along_vector(components, nbins=10, vector=[1,1]):
    """
    project data on line defined by unit vector and bin it
    """
    ndim  = len(vector)
    u_vec = vector/np.linalg.norm(vector)
    projection = np.dot(components[:,0:ndim],u_vec)
    data_percentile = np.array([stats.percentileofscore(projection, a) for a in projection])
    binsize = int(100/nbins)
    bins_percentile = np.arange(0, 100+binsize, binsize)
    data_binned_indices = np.digitize(data_percentile, bins_percentile, right=True)
    index = assignment_to_index(data_binned_indices, nbins+1, istart=1)
    return data_binned_indices, index

def bin_component(components, nbins=10, i_component=0):
    """
    bin along i_component
    """
    data_percentile = np.array([stats.percentileofscore(components[:,i_component], a) for a in components[:,i_component]])
    binsize = int(100/nbins)
    bins_percentile = np.arange(0, 100+binsize, binsize)
    data_binned_indices = np.digitize(data_percentile, bins_percentile, right=True)
    index = assignment_to_index(data_binned_indices, nbins+1, istart=1)
    return data_binned_indices, index

## GMM

def scan_gmm(V,n_components=10,plot=True, do_return=False):
    """
    """
    n_components_range = range(1, n_components+1)
    lowest_bic = np.infty
    bic = []
    for i_components in n_components_range:
        gmm = mixture.GaussianMixture(n_components=i_components)
        gmm.fit(V)
        bic.append(gmm.bic(V))
        print(" BIC(",i_components,") = ",bic[-1],end="")
        if bic[-1] < lowest_bic:
            lowest_bic = bic[-1]
            best_gmm = gmm
    bic = np.array(bic)
    if plot:
        plot_vector(bic) 
    else:
        print('BIC')
        print(bic)
    if do_return:
        return best_gmm,bic

def assign_gmm(V,n_components=2):
    """
    """
    gmm = mixture.GaussianMixture(n_components=n_components)
    gmm.fit(V)
    assignment = gmm.predict(V)
    index = assignment_to_index(assignment, n_components)
    return assignment,index

## SEGMENTATION

def segment_map(input_map, length_scale=1):
    """
    """
    distance = ndimage.distance_transform_edt(input_map)
    local_maxi = peak_local_max(distance, 
                                indices=False, 
                                footprint=np.ones((length_scale,length_scale,length_scale)), 
                                labels=input_map)
    markers = morphology.label(local_maxi)
    labelled_map = watershed(-distance, markers, mask=input_map)
    return labelled_map

## TOOLS

def assignment_to_index(assignment, nbins, istart=0):
    index = []
    for i in np.arange(istart,nbins):
        index.append(np.ma.where(assignment == i))
        print(". component ",i," has ",index[-1][0].shape[0]," particles")
    return index

## PLOTS

def plot_vector(x):
    """
    """
    fig = plt.figure(figsize=(2,2))
    plt.plot(np.arange(1,x.shape[0]+1),x,marker='o',color='black')

#
