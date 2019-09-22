import numpy as np
import re
from matplotlib import pyplot as plt
from scipy import linalg
from sklearn.preprocessing import normalize
#########################
# interface with RELION #
#########################
def get_cryosparc_projections(proj_file, nmodes=3, norm=10000):
    """
    """
    print('Load coordinates:')
    data = np.array(np.load(proj_file).tolist())
    i = 2*nmodes-1
    components = data[:,-i::2]/norm
    print('Number of components :',components.shape[1])
    print('Number of particles  :',components.shape[0])
    return components

def get_relion_projections(proj_file,eigvec_file, as_is=False):
    """
    get_relion_projections
    """
    print('Load data [check out eigenvectors and eigenvalues]:')
    ids,dims,U,L,Vt = load_pca(proj_file,eigvec_file, as_is=as_is)
    return Vt.T

def split_relion_projections(proj_file, eigvec_file, c=None, keyword='split_'):
    """
    """
    ids,dims,U,L,Vt = load_pca(proj_file,eigvec_file)
    save_cluster(ids,dims,U,L,Vt,c=c,keyword=keyword)

def load_pca(proj_file,eigvec_file,as_is=True,plot=True):
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
    pfile = np.genfromtxt(proj_file,dtype='str')
    ids   = pfile[:,0]
    proj, evalue  = normalize(pfile[:,1:],axis=0,return_norm=True)
    # load eigenvectors
    vfile = np.genfromtxt(eigvec_file,dtype='str')
    dims = vfile[0,:]
    evector = vfile[1:,:].astype(np.float)
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
            plot_eig(dims,U.T,L)
        return ids,dims,U,L,Vt
    else:
        print('ERROR in dimensions...')

def save_cluster(ids,dims,U,L,Vt,c=None,keyword='save_'):
    """
    """
    if c is not None:
        # save eigenvector
        filename=keyword+'_eigvec.dat'
        f = open(filename,'w')
        for item in dims:
            f.write("%s " % item)
        f.write("\n")
        for i_component in np.arange(0,U.shape[1]):
            comp = U[:,i_component]
            for item in comp:
                f.write("%e " % item)
            f.write("\n")
        f.close()
        # for each cluster, save projections, scaled by eigenvalue
        projs = np.dot(np.diag(L),Vt)
        n_components = np.max(np.unique(c))+1
        for i in np.arange(0,n_components):
            filename=keyword+'_cluster_'+str(i)+'_projs.dat'
            f = open(filename,'w')
            index = np.ma.where(c == i)
            ids_kept = ids[index]
            prj_kept = projs[:,index[0]]
            for j in np.arange(0,len(index[0])):
                f.write("%s " % ids_kept[j])
                prjlist = prj_kept[:,j]
                for item in prjlist:
                    f.write("%f " % item)
                f.write("\n")
            f.close()

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

###
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

###
def plot_eig(dims,eigvec,eigval,figsize=4):
    """
    """
    dimarray = np.arange(len(dims))
    fig, (ax1,ax2) = plt.subplots(1,2,sharey=True,figsize=(2*figsize,figsize))
    im1 = ax1.imshow(eigvec.T,cmap='seismic')
    ax1.set_ylabel('component #')
    ax1.set_yticks(dimarray)
    ax1.set_yticklabels(np.arange(1,len(dims)+1))
    ax1.set_xlabel('bodies parameters')
    ax1.set_xticks(dimarray)
    ax1.set_xticklabels(dims,rotation='vertical')
    cbar = ax1.figure.colorbar(im1, ax=ax1)
    im2 = ax2.barh(dimarray,eigval,color='black')
    ax2.set_yticks(dimarray)
    ax2.set_yticklabels(np.arange(1,len(dims)+1))
    ax2.set_xlabel('eigenvalue')
    plt.tight_layout()
    plt.show()
