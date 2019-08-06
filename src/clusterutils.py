import numpy as np
from matplotlib import pyplot as plt
from sklearn import mixture

def scan_gmm(V,n_components=3,plot=True):
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
    return best_gmm,bic

def assign_gmm(V,n_components=2):
    """
    """
    gmm = mixture.GaussianMixture(n_components=n_components)
    gmm.fit(V)
    assignment = gmm.predict(V)
    index = []
    for i in np.arange(0,n_components):
        index.append(np.ma.where(assignment == i))
        print(". component ",i," has ",index[-1][0].shape[0]," particles")
    return assignment,index

#########
# < simple plots

def plot_vector(x):
    """
    """
    fig = plt.figure(figsize=(2,2))
    plt.plot(np.arange(1,x.shape[0]+1),x,marker='o',color='black')

#
