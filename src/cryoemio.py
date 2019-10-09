import numpy as np
import h5py
import os
import mrcfile
import CifFile
from CifFile import ReadCif

## RELION related ##

def star_reader(filename):
    """
    star_reader: reads a star file and returns a dictionary
    """
    key='data_'
    new_key='data_metadata'
    new_starfile='.tmp.star'
    #
    with open(filename) as f1:
        with open(new_starfile, 'w') as f2:
            lines = f1.readlines()
            for line in lines:
                if(line.startswith(key)):
                    f2.write(new_key)
                else:
                    f2.write(line)
    data = ReadCif(new_starfile, grammar='STAR2')
    os.remove(new_starfile)
    print('Number of particles in this star file: {}'.format(len(data['metadata'][data['metadata'].keys()[0]])))
    print("The entries in the returned dictionary are:")
    print("data['metadata'].keys(): {}".format(data['metadata'].keys()))
    return data

def star2np(filename, key, keytype=float):
    """
    star2np: returns an item from a star file as a numpy array
    """
    data  = star_reader(filename)
    array = starkey2np(data, key, keytype=keytype)
    return array

def starkey2np(data, key, keytype=float):
    """
    starkey2np: star file has already been read.
    """
    if key in data['metadata']:
        array = np.array(data['metadata'][key]).astype(keytype)
    return array

def star2hdf5(starfile, hdf5file, path_to_mrcs):
    print("> Reading star file...")
    data = star_reader(starfile)
    n_particles = len(data['metadata']['_rlnimagename'])
    image = rlnimagename2image(data['metadata']['_rlnimagename'][0], path_to_mrcs)
    nx,ny = image.shape
    image_stack = np.empty((n_particles, nx, ny))
    print("> Reading corresponding particle mrcs files...")
    for i in np.arange(n_particles):
        image_stack[i,...] = rlnimagename2image(data['metadata']['_rlnimagename'][i], path_to_mrcs)
    with h5py.File(hdf5file, 'w') as hf:
        print("> Writing particle images to h5 file...")
        hf.create_dataset('particles', data=image_stack)
        print("> Writin metadata to h5 file...")
        for key in data['metadata'].keys():
            try:
                content = np.array(data['metadata'][key], dtype='float_')
            except:
                content = np.array(data['metadata'][key], dtype='<S')
            hf.create_dataset(key, data=content)
    print('> Done! Output at {}'.format(hdf5file))

def star2hdf5_serial(starfile, hdf5file, path_to_mrcs, istep=-1):
    print("> Reading star file...")
    data = star_reader(starfile)
    n_particles = len(data['metadata']['_rlnimagename'])
    if(istep<=0):
        istep = int(n_particles/100)
    with h5py.File(hdf5file, 'w') as hf:
        print("> Writin metadata to h5 file...")
        for key in data['metadata'].keys():
            print('>>> {}'.format(key))
            try:
                content = np.array(data['metadata'][key], dtype='float_')
            except:
                content = np.array(data['metadata'][key], dtype='<S')
            hf.create_dataset(key, data=content)
        print('> Writing {} particle images to h5 file...'.format(n_particles))
        for i in np.arange(n_particles):
            if(i % 10 == 0):
                print('[{}/{}] '.format(i, n_particles), end=" ")
            image = rlnimagename2image(data['metadata']['_rlnimagename'][i], path_to_mrcs)
            image = image[np.newaxis,...]
            if(i==0):
                hf.create_dataset('particles', data=image, maxshape=(None,image.shape[1],image.shape[2]), chunks=True)
            else:
                hf['particles'].resize( (hf['particles'].shape[0] + image.shape[0]), axis=0 )
                hf['particles'][-image.shape[0]:] = image
        print('')
    print('> Done! Output at {}'.format(hdf5file))

def rlnimagename2image(particle, path_to_mrcs):
    frame, relpath = particle.split('@')
    iframe = np.int(frame)
    mrcs_filename = relpath.split('/')[-1]
    mrcs_file = path_to_mrcs+mrcs_filename
    mrcs_data = mrc2data(mrc_file = mrcs_file)
    image = mrcs_data[iframe-1,...]
    return image

def get_mrcorder_list(data, path):
    """
    the order is different when:
    - read from the star file (star2hdf5_serial)
    - read from the list glob (mrclist2hdf5)
    this function provides the index order
    """
    list_in_data = []
    for image in data['_rlnimagename']:
        string = str(image).split("'")[1]
        frame, relpath = string.split('@')
        fname  = relpath.split('/')[-1]
        list_in_data.append(fname)
    #
    list_in_path = []
    for mrcs in glob.glob(path+'*.mrcs'):
        fname = mrcs.split('/')[-1]
        list_in_path.append(fname)
    #
    index_order = []
    for key in list_in_path:
        hits = [i for i, j in enumerate(list_in_data) if j == key]
        for hit in hits:
            index_order.append(hit)
    #
    return index_order

##

def simio(pdbdir, pdb_keyword, mrcdir, mrc_keyword, action='define'):
    """ simio
    """
    pdb_file = pdbdir+pdb_keyword+'.pdb'
    crd_file = mrcdir+mrc_keyword+'.txt'
    mrc_file = mrcdir+mrc_keyword+'.mrc'
    log_file = mrcdir+mrc_keyword+'.log'
    inp_file = mrcdir+mrc_keyword+'.inp'
    h5_file  = mrcdir+mrc_keyword+'.h5'
    if(action=='define'):
        return pdb_file, mrc_file, crd_file, log_file, inp_file, h5_file
    elif(action=='clean'):
        for f in (mrc_file, crd_file, log_file, inp_file):
            if os.path.isfile(f):
                os.remove(f)

def mrc2data(mrc_file = None):
    """ mrc2data
    """
    if mrc_file is not None:
        with mrcfile.open(mrc_file, 'r+', permissive=True) as mrc:
            micrograph = mrc.data
        if(len(micrograph.shape)==2):
            micrograph = micrograph[np.newaxis,...]
        return micrograph

###### / dictionary approach
def mrc2dic2hdf5(mrc_file = None, h5_file = None, dic = None):
    """ mrc2hdf5
    """
    if mrc_file is not None:
        micrograph = mrc2data(mrc_file = mrc_file)
        if dic is None:
            dic = {}
        dic['data'] = micrograph
        save_dict_to_hdf5(dic, h5_file)

def data_and_dic_2hdf5(data, h5_file, dic = None):
    """
    """
    if dic is None:
        dic = {}
    dic['data'] = data
    save_dict_to_hdf5(dic, h5_file)

############# / specific to TEM-simulator
def add_crd_to_dic(crd_file = None, dic = None):
    """ add_crd_to_dic: [specific to TEM-simulator input/output]
    this function reads the text content of a coordinates file and adds it to the input dictionary.
    Returns the richer dictionary.
    """
    if dic is None:
        dic = {}
    if crd_file is not None:
        crd = np.genfromtxt(crd_file, skip_header=3)
        dic['coordinates'] = crd
    return dic
#
def add_crd_to_h5(input_h5file = None, input_crdfile = None, output_h5file = None):
    """ add_crd_to_h5: [specific to TEM-simulator input/output]
    this function assumes an .hdf5 file containing particle images in a dictionary. Potentially some other info.
    We want to add a field here where the particle coordinates read from the crd file are added.
    """
    if input_h5file is not None:
        dic = load_dict_from_hdf5(input_h5file)
    else:
        dic = {}
    if input_crdfile is not None:
        crd = np.genfromtxt(input_crdfile, skip_header=3)
        dic['coordinates'] = crd
        save_dict_to_hdf5(dic, output_h5file)
    else:
        print("Error")
############ specific to TEM-simulator / #####
##### dictionary approach / ####

def mrclist2hdf5( mrc_list=None, h5_file=None , verbose=False):
    """
    """
    if h5_file is not None:
        with h5py.File(h5_file, 'w') as hf:
            i=0
            if mrc_list is not None:
                for mrc in mrc_list:
                    new_data = mrc2data(mrc_file = mrc)
                    if verbose:
                        print('{0} {1}'.format(mrc,new_data.shape))
                    if(i==0):
                        hf.create_dataset('particles', data=new_data,maxshape=(None,new_data.shape[1],new_data.shape[2]), chunks=True)
                    else:
                        hf['particles'].resize( (hf['particles'].shape[0] + new_data.shape[0]), axis=0 )
                        hf['particles'][-new_data.shape[0]:] = new_data
                    i+=1

def loadhdf5(h5_file, verbose=False):
    """ 
    ...
    """
    data = h5py.File(h5_file,'r')
    if verbose:
        print('List of datasets in .hdf5 file:')
        data.visititems(print_attrs)
    return data

def print_attrs(name, obj):
    print( name )
    for key, val in obj.attrs.items():
        print("    %s: %s" % (key, val))


# The following hdf5 save and load functions are taken from here: 
#https://codereview.stackexchange.com/questions/120802/recursively-save-python-dictionaries-to-hdf5-files-using-h5py
def save_dict_to_hdf5(dic, filename):
    """
    ....
    """
    with h5py.File(filename, 'w') as h5file:
        recursively_save_dict_contents_to_group(h5file, '/', dic)
#
def recursively_save_dict_contents_to_group(h5file, path, dic):
    """
    ....
    """
    for key, item in dic.items():
        if isinstance(item, (np.ndarray, np.int64, np.float64, int, float, str, bytes)):
            h5file[path + key] = item
        elif isinstance(item, type(None)):
            h5file[path + key] = str('None')
        elif isinstance(item, dict):
            recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
        else:
            raise ValueError('Cannot save %s type'%type(item))
#
def load_dict_from_hdf5(filename):
    """
    ....
    """
    with h5py.File(filename, 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')
#
def recursively_load_dict_contents_from_group(h5file, path):
    """
    ....
    """
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item[()] #item.value
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans
####################
