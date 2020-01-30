import numpy as np
import re
from matplotlib import pyplot as plt
from scipy import linalg
from sklearn.preprocessing import normalize

## INTERFACES

def get_cryosparc_projections(proj_file, nmodes=3, norm=10000):
    """
    get_cryosparc_projections: read passthrough file, outputs components array
    +++++++++++++++++++++++++
    """
    print('Load coordinates:')
    data = np.array(np.load(proj_file).tolist())
    i = 2*nmodes-1
    components = data[:,-i::2]/norm
    print('Number of components :',components.shape[1])
    print('Number of particles  :',components.shape[0])
    return components

def get_relion_projections(proj_file,eigvec_file, as_is=False, figname=''):
    """
    get_relion_projections: read output from multibody refinement, outputs components array
    ++++++++++++++++++++++
    """
    print('Load data [check out eigenvectors and eigenvalues]:')
    ids,dims,U,L,Vt = load_pca(proj_file,eigvec_file, as_is=as_is, figname=figname)
    return Vt.T

def split_relion_projections(proj_file, eigvec_file, c=None, keyword='split_'):
    """
    split_relion_projections: read output from multibody refinement and given assignment,
    ++++++++++++++++++++++++  splits dataset accordingly and writes each cluster out.
    """
    ids,dims,U,L,Vt = load_pca(proj_file,eigvec_file, plot=False)
    save_cluster(ids,dims,U,L,Vt,c=c,keyword=keyword)

## PCA

def load_pca(proj_file,eigvec_file,as_is=True,plot=True, figname=''):
    """ load_pca:
    INPUTS
    ------
        . proj_file  : output file from flex_analyze containing particle coordinates in PC space
        . eigvec_file: output file from flex_analyze containing the definition of the eigen vectors
        . as_is      : if False, reconstitute the data first, and output SVD. 
                       otherwise, consider that input is OK
    RETURNS
    -------
        . ids[nsample]         : particle names or IDs
        . dims[n_dim]          : names of the body dimensions (trans-rot) 
        . evector[n_dim,n_cmp] : columns are eigenvectors 
        . evalue[n_cmp]        : associated eigenvalues
        . proj[n_cmp,n_sample] : each column is the coordinate of a particle
    """
    # load projections and extract eigenvalues
    proj, evalue, ids = read_relion_projection_file(proj_file)
    # load eigenvectors
    evector, dims = read_relion_eigenvector_file(eigvec_file)
    # check dimensions
    if proj.shape[1] == evector.shape[0] == evalue.shape[0] :
        if as_is:
            U  = evector.T
            L  = evalue
            Vt = proj.T
        else:
            data = recompose(evector.T,evalue,proj.T)
            demeaned_data = de_mean(data)
            U,L,Vt = get_svd(demeaned_data.T)
        print('Number of dimensions :',dims.shape[0])
        print('Number of components :',Vt.shape[0])
        print('Number of particles  :',Vt.shape[1])
        if plot :
            plot_eig(dims,U.T,L,figname=figname)
        return ids,dims,U,L,Vt
    else:
        print('ERROR in dimensions...')

def get_svd(data):
    """
    """
    Vsvd,Lsvd,Usvd = linalg.svd(data,full_matrices=False)
    return Usvd.T,Lsvd,Vsvd.T

def recompose(U,L,Vt):
    """ recompose
    """
    data = np.dot(U,np.dot(np.diag(L),Vt))
    return data

def de_mean(data):
    """
    """
    return (data.T - np.mean(data, axis=1)).T


## ICA stuff

def ave_score(X,n,niter=100,fun='logcosh'):
    # notice that X components are column-oriented
    score_ave=[]
    score_var=[]
    for i in np.arange(0,n,1):
        score_tmp = []
        for j in np.arange(0,niter,1):
            score_tmp.append(negent_score(X[:,i],fun))
        score_ave.append(np.mean(score_tmp))
        score_var.append(np.var(score_tmp))
    return score_ave,score_var

def negent_score(X,fun='logcosh'):
    # We compute J(X) = [E(G(X)) - E(G(Xgauss))]**2
    # We consider X (and Xgauss) to be white, in the sense that E(X,X.T)=I
    # The expectation being approximated by the sample mean in our case: np.dot(X,X.T)/n=I
    # In practice, we assume that X has already been normalized by its length [np.dot(X,X.T)=I]
    # so we rescale by np.sqrt(n) before we take the expectation value of G(X).
    length=len(X)
    Xscale = X*np.sqrt(length)
    Xgauss = np.random.randn(length)
    if(fun == 'logcosh'):
        n1 = np.mean(f_logcosh(Xscale)) #np.sum(f_logcosh(Xscale))
        n2 = np.mean(f_logcosh(Xgauss)) #np.sum(f_logcosh(Xgauss))
    elif(fun == 'exp'):
        n1 = np.mean(f_exp(Xscale))     #np.sum(f_exp(Xscale))
        n2 = np.mean(f_exp(Xgauss))     #np.sum(f_exp(Xgauss))
    elif(fun == 'rand'):
        n1 = np.mean(f_logcosh(Xgauss)) #np.sum(f_logcosh(Xgauss))
        n2 = 0
    negent = (n2-n1)**2
    return negent

def f_logcosh(X):
    return np.log(np.cosh(X))

def f_exp(X):
    return -np.exp(-(X**2)/2)

## READ-WRITE

def save_cluster(ids,dims,U,L,Vt,c=None,keyword='save_'):
    """
    """
    if c is not None:
        # save eigenvector
        filename=keyword+'_eigvec.dat'
        print('.. writing {0}'.format(filename))
        write_relion_eigenvector_file(filename, U, dims)
        # for each cluster, save projections, scaled by eigenvalue
        projs = np.dot(np.diag(L),Vt)
        n_components = np.max(np.unique(c))+1
        for i in np.arange(0,n_components):
            index = np.ma.where(c == i)
            filename=keyword+'_cluster_'+str(i)+'_projs.dat'
            print('.. writing {0}'.format(filename))
            write_relion_projection_file(filename, projs, ids, index)

def read_relion_projection_file(proj_file):
    pfile = np.genfromtxt(proj_file,dtype='str')
    ids   = pfile[:,0]
    proj, evalue  = normalize(pfile[:,1:],axis=0,return_norm=True)
    return proj, evalue, ids

def write_relion_projection_file(filename, projs, ids, index):
    f = open(filename,'w')
    ids_kept = ids[index]
    prj_kept = projs[:,index[0]]
    for j in np.arange(0,len(index[0])):
        f.write("%s " % ids_kept[j])
        prjlist = prj_kept[:,j]
        for item in prjlist:
            f.write("%f " % item)
        f.write("\n")
    f.close()

def read_relion_eigenvector_file(eigvec_file):
    vfile = np.genfromtxt(eigvec_file,dtype='str')
    dims = vfile[0,:]
    evector = vfile[1:,:].astype(np.float)
    return evector, dims

def write_relion_eigenvector_file(filename, evector, dims):
    f = open(filename,'w')
    for item in dims:
        f.write("%s " % item)
    f.write("\n")
    for i_component in np.arange(0,evector.shape[1]):
        comp = evector[:,i_component]
        for item in comp:
            f.write("%e " % item)
        f.write("\n")
    f.close()

### STAR 

def dat2star(starfile,datfile,outfile='test.star'):
    particles_ids  = particles_read(starfile)
    cluster_ids    = particles_read(datfile,filetype='dat')
    particles_mask = np.in1d(particles_ids, cluster_ids, assume_unique=True)
    particles_write(starfile,mask=particles_mask,output_file=outfile)

def particles_read(file,filetype='star'):
    icol=5
    if filetype == 'dat':
        icol = 0
    data = open(file,"r")
    ids = []
    for line in data:
        if re.search('@', line) :
            ids.append(line.split()[icol])
    print('Total number of particles = ',len(ids))
    return np.asarray(ids)

def particles_write(file,mask=None,output_file='test.dat'):
    data   = open(file,"r")
    output = open(output_file,"w")
    n_particles=0
    for line in data:
        if re.search('@', line) :
            n_particles+=1
            if mask[n_particles-1]:
                output.write(line)
        elif re.search('data_images_body_', line):
            break
        else:
            output.write(line)
            n_particles=0
    output.close()
    print('Wrote particles to file ',output_file)

## PLOT

def plot_eig(dims,eigvec,eigval,figsize=4, show_explained_variance=True,figname=''):
    """
    """
    dimarray = np.arange(len(dims))
    #
    fig, (ax1,ax2) = plt.subplots(1,2,sharey=True,figsize=(2*figsize,figsize))
    #
    im1 = ax1.imshow(eigvec.T,vmin=-1,vmax=1,cmap='seismic')
    ax1.set_ylabel('component #')
    ax1.set_yticks(dimarray)
    ax1.set_yticklabels(np.arange(1,len(dims)+1))
    ax1.set_xlabel('bodies parameters')
    ax1.set_xticks(dimarray)
    ax1.set_xticklabels(dims,rotation='vertical')
    cbar = ax1.figure.colorbar(im1, ax=ax1)
    #
    value = eigval
    if show_explained_variance:
        value = (eigval**2)/np.sum(eigval**2)
    dimarray_display = np.arange(-1,len(dims)+1)
    eigval_display = np.append(np.insert(value,0,0),0)
    im2 = ax2.barh(dimarray_display,eigval_display,height=0.8, color='black')
    #ax2.set_yticks(dimarray_display)
    #ax2.set_yticklabels(np.arange(0,len(dims)+2))
    if show_explained_variance:
        ax2.set_xlabel('ratio of explained variance')
    else:
        ax2.set_xlabel('eigenvalue')
    plt.tight_layout()
    plt.show()
    if(figname):
        fig.savefig(figname)
