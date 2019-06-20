import os
import sys
import subprocess
import numpy as np
import mdtraj as md
from matplotlib import pyplot as plt
import math
from scipy.stats import special_ortho_group
from scipy.spatial.distance import pdist
#
def isRotationMatrix(R) :
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6

def rotationMatrixToEulerAngles(R) :
    assert(isRotationMatrix(R))
    sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
    singular = sy < 1e-6
    if  not singular :
        x = math.atan2(R[2,1] , R[2,2])
        y = math.atan2(-R[2,0], sy)
        z = math.atan2(R[1,0], R[0,0])
    else :
        x = math.atan2(-R[1,2], R[1,1])
        y = math.atan2(-R[2,0], sy)
        z = 0
    x = (x*180)/np.pi
    y = (y*180)/np.pi
    z = (z*180)/np.pi
    return np.array([x, y, z])

def get_rotlist(numpart):
    """ get_rotlist
    """
    rotlist = []
    for x in range(0,numpart+1):
        x = special_ortho_group.rvs(3)
        y = rotationMatrixToEulerAngles(x)
        rotlist.append(y)
    return rotlist
        
def write_crd_file(numpart, xrange=np.arange(-100,110,10), yrange=np.arange(-100,110,10), crd_file='crd.txt'):
    """ write_crd_file
    The table should have 6 columns. The first three columns are x, y, and z coordinates of
    the particle center. The following three columns are Euler angles for rotation around
    the z axis, then around the x axis, and again around the z axis. 
    Coordinates are in nm units, and angles in degrees.
    """
    if os.path.exists(crd_file):
         print(crd_file+" already exists.")          
    else:
        rotlist  = get_rotlist(numpart)
        crd = open(crd_file, "w")
        crd.write('# File created by TEM-simulator, version 1.3.\n')
        crd.write('{numpart}  6\n'.format(numpart=numpart))
        crd.write('#            x             y             z           phi         theta           psi  \n')
        l = 0
        for y in yrange:
            for x in xrange:
                if l == int(numpart):
                    break
                crd_table = {'x': x, 'y': y, 'z': 0, 'phi': rotlist[l][0], 'theta': rotlist[l][1], 'psi': rotlist[l][2]}
                crd.write('{0[x]:14.4f}{0[y]:14.4f}{0[z]:14.4f}{0[phi]:14.4f}{0[theta]:14.4f}{0[psi]:14.4f}\n'.format(crd_table))
                l += 1
        crd.close()
#
def fill_parameters_dictionary(mrc_file = None,
                               pdb_file = None, voxel_size=0.1, particle_name='toto', particle_mrcout=None,
                               crd_file = None, sample_dimensions=[1200,50,150],
                               beam_params=[300,1.3,100,0], dose=None,
                               optics_params=[81000,2.7,2.7,50,3.5,0.1,1.0,0,0], defocus=None, optics_defocout=None,
                               detector_params=[5760,4092,5,32,'yes',0.5,0.0,0.0,1.0,0,0], noise=None,
                               log_file='simulator.log', seed = -1234):
    """ fill_parameters_dictionary
        
        Parameters: default values in parentheses
        -----------
        - mrc_file (output) [MANDATORY]: Micrograph file
        *** PARTICLE ***
        - pdb_file (input) [MANDATORY]: PDB file of sample
        - voxel size (0.1)      : The size of voxels in the particle map in nm.
        - particle_name ('toto'): Name of the particle. Not very important.
        - particle_mrcout (None): if not None, volume map of sample is written.
        *** GRID ***
        - crd_file (input) [MANDATORY]: Coordinates of the sample copies
        - sample_dimensions:
            . index 0 (1200): diameter in nm
            . index 1 (50)  : thickness at center in nm
            . index 2 (150) : thickness at edge in nm.
        *** MICROSCOPE ***
        - beam_params:
            . index 0 (300)    : voltage in kV
            . index 1 (1.3)    : energy spread in V
            . index 2 (100): dose per image in e/nm**2
            . index 3 (0)      : standard deviation of dose per image
        - dose (None): if not None, overrides beam_params[2]
        - optics_params:
            . index 0 (81000) : magnification (81000; 105000; 130000)
            . index 1 (2.7)   : spherical aberration in mm
            . index 2 (2.7)   : chromatic aberration in mm
            . index 3 (50)    : diameter in um of aperture in back focal plane (50-100)
            . index 4 (3.5)   : focal length in mm of primary lens
            . index 5 (0.1)   : aperture angle in mrad of the beam furnished by 
                                the condenser lens
            . index 6 (1.0)   : nominal defocus value in um
            . index 7 (0)     : standard deviation of a systematic error added 
                                to the nominal defocus, measured in um. 
                                The same error is added to the defocus of every image.
            . index 8 (0)     : standard deviation of a nonsystematic error added 
                                to the nominal defocus and the systematic error, 
                                measured in um. 
                                A new value of the error is computed for every image
        - defocus (None): if not None, overrides optics_params[6]
        - optics_defocout (None): if not None, defocus values written to file
        *** DETECTOR ***
        - detector_params: 
            . index 0  (5760) : number of pixels on detector along x axis
            . index 1  (4092) : number of pixels on detector along y axis
            . index 2  (5)   : physical pixel size in um
            . index 3  (32)   : detector gain: average number of counts per electron
            . index 4  ('yes'): quantized electron waves result in noise
            . index 5  (0.5)  : detector quantum efficiency
            . index 6  (0)  : parameter of MTF
            . index 7  (0)  : parameter of MTF
            . index 8  (1)  : parameter of MTF
            . index 9  (0)   : parameter of MTF
            . index 10 (0)   : parameter of MTF
        - noise (None): if not None, overrides detector_params[4]
        *** MISC ***
        - log_file ('simulator.log'): Log file for the run
        - seed (-1234): seed for the run       
    """
    # see if we need overrides
    if dose is not None:
        beam_params[2]     = dose
    if defocus is not None:
        optics_params[6]   = defocus
    if noise is not None:
        detector_params[4] = noise
    # fill the dictionary
    dic = {}
    dic['simulation'] = {}
    dic['simulation']['seed']    = seed
    dic['simulation']['logfile'] = log_file
    dic['sample'] = {}
    dic['sample']['diameter']         = sample_dimensions[0] # diameter in nm
    dic['sample']['thickness_center'] = sample_dimensions[1] # thickness at center in nm
    dic['sample']['thickness_edge']   = sample_dimensions[2] # thickness at edge in nm
    dic['particle'] = {}
    dic['particle']['name']       = particle_name
    dic['particle']['voxel_size'] = voxel_size
    dic['particle']['pdb_file']   = pdb_file
    if particle_mrcout is None:
        dic['particle']['map_file_re_out'] = None
    else:
        key = mrc_file.split('.mrc')[0]
        dic['particle']['map_file_re_out'] = key+'_real.mrc'
        dic['particle']['map_file_im_out'] = key+'_imag.mrc'
    dic['particleset'] = {}
    dic['particleset']['name']     = particle_name
    dic['particleset']['crd_file'] = crd_file
    dic['beam'] = {}
    dic['beam']['voltage']      = beam_params[0] # voltage in kV
    dic['beam']['spread']       = beam_params[1] # energy spread in V
    dic['beam']['dose_per_im'] = beam_params[2] # dose per image in e/nm**2
    dic['beam']['dose_sd']      = beam_params[3] # standard deviation of dose per image
    dic['optics'] = {}
    dic['optics']['magnification']         = optics_params[0] # magnification
    dic['optics']['cs']                    = optics_params[1] # spherical aberration in mm
    dic['optics']['cc']                    = optics_params[2] # chromatic aberration in mm
    dic['optics']['aperture']              = optics_params[3] # diameter in um of aperture in back focal plane
    dic['optics']['focal_length']          = optics_params[4] # focal length in mm of primary lens
    dic['optics']['cond_ap_angle']         = optics_params[5] # aperture angle in mrad of the beam furnished by the condenser lens
    dic['optics']['defocus_nominal']       = optics_params[6] # nominal defocus value in um
    dic['optics']['defocus_syst_error']    = optics_params[7] # standard deviation of a systematic error added to the nominal defocus, measured in um. The same error is added to the defocus of every image
    dic['optics']['defocus_nonsyst_error'] = optics_params[8] # standard deviation of a nonsystematic error added to the nominal defocus and the systematic error, measured in um. A new value of the error is computed for every image
    if optics_defocout is None:
        dic['optics']['defocus_file_out']  = None
    else:
        dic['optics']['defocus_file_out']  = optics_defocout  # file to which defocus values are written
    dic['detector'] = {}
    dic['detector']['det_pix_x']        = detector_params[0] # number of pixels on detector along x axis
    dic['detector']['det_pix_y']        = detector_params[1] # number of pixels on detector along y axis
    dic['detector']['pixel_size']       = detector_params[2] # physical pixel size in um
    dic['detector']['gain']             = detector_params[3] # detector gain: average number of counts per electron
    dic['detector']['use_quantization'] = detector_params[4] # quantized electron waves result in noise
    dic['detector']['dqe']              = detector_params[5] # detector quantum efficiency
    dic['detector']['mtf_a']            = detector_params[6] # parameter of MTF
    dic['detector']['mtf_b']            = detector_params[7] # parameter of MTF
    dic['detector']['mtf_c']            = detector_params[8] # parameter of MTF
    dic['detector']['mtf_alpha']        = detector_params[9] # parameter of MTF
    dic['detector']['mtf_beta']         = detector_params[10] #parameter of MTF
    dic['detector']['image_file_out']   = mrc_file # file with resulting micrograph
    return dic               
                
def write_inp_file(inp_file='input.txt', dict_params=None):
    """ write_inp_file
    """
    if dict_params is not None:
        inp = open(inp_file, "w")
        inp.write('=== simulation ===\n'
              'generate_micrographs = yes\n'
              'rand_seed = {0[seed]}\n'
              'log_file = {0[logfile]}\n'.format(dict_params['simulation']))
        inp.write('=== sample ===\n'
              'diameter = {0[diameter]:d}\n'
              'thickness_edge = {0[thickness_edge]:d}\n'
              'thickness_center = {0[thickness_center]:d}\n'.format(dict_params['sample']))
        inp.write('=== particle {0[name]} ===\n'
              'source = pdb\n'
              'voxel_size = {0[voxel_size]}\n'
              'pdb_file_in = {0[pdb_file]}\n'.format(dict_params['particle']))
        if dict_params['particle']['map_file_re_out'] is not None:
            inp.write('map_file_re_out = {0[map_file_re_out]}\n'
                  'map_file_im_out = {0[map_file_im_out]}\n'.format(dict_params['particle']))
        inp.write('=== particleset ===\n'
              'particle_type = {0[name]}\n'
              'particle_coords = file\n'
              'coord_file_in = {0[crd_file]}\n'.format(dict_params['particleset']))
        inp.write('=== geometry ===\n'
              'gen_tilt_data = yes\n'
              'tilt_axis = 0\n'
              'ntilts = 1\n'
              'theta_start = 0\n'
              'theta_incr = 0\n'
              'geom_errors = none\n')
        inp.write('=== electronbeam ===\n'
              'acc_voltage = {0[voltage]}\n'
              'energy_spread = {0[spread]}\n'
              'gen_dose = yes\n'
              'dose_per_im = {0[dose_per_im]}\n'
              'dose_sd = {0[dose_sd]}\n'.format(dict_params['beam']))
        inp.write('=== optics ===\n'
              'magnification = {0[magnification]}\n'
              'cs = {0[cs]}\n'
              'cc = {0[cc]}\n'
              'aperture = {0[aperture]}\n'
              'focal_length = {0[focal_length]}\n'
              'cond_ap_angle = {0[cond_ap_angle]}\n'
              'gen_defocus = yes\n'
              'defocus_nominal = {0[defocus_nominal]}\n'
              'defocus_syst_error = {0[defocus_syst_error]}\n'
              'defocus_syst_error = {0[defocus_nonsyst_error]}\n'.format(dict_params['optics']))
        if dict_params['optics']['defocus_file_out'] is not None:
            inp.write('defocus_file_out = {0[defocus_file_out]}\n'.format(dict_params['optics']))
        inp.write('=== detector ===\n'
              'det_pix_x = {0[det_pix_x]}\n'
              'det_pix_y = {0[det_pix_y]}\n'
              'pixel_size = {0[pixel_size]}\n'
              'gain = {0[gain]}\n'
              'use_quantization = {0[use_quantization]}\n'
              'dqe = {0[dqe]}\n'
              'mtf_a = {0[mtf_a]}\n'
              'mtf_b = {0[mtf_b]}\n'
              'mtf_c = {0[mtf_c]}\n'
              'mtf_alpha = {0[mtf_alpha]}\n'
              'mtf_beta = {0[mtf_beta]}\n'
              'image_file_out = {0[image_file_out]}\n'.format(dict_params['detector']))
        inp.close()
#
def get_xyz_from_pdb(filename=None):
    """ get_xyz_from_pdb
    """
    if filename is not None:
        traj = md.load(filename)
        atom_indices = traj.topology.select('name CA or name P')
        traj_small = traj.atom_slice(atom_indices)
        return traj_small.xyz

def get_Dmax(filename=None):
    """ get_Dmax
    """
    if filename is not None:
        xyz = get_xyz_from_pdb(filename)
        distance= pdist(xyz[0,...])
        return np.amax(distance)
#
def define_grid_in_fov(sample_dimensions, optics_params, detector_params, pdb_file=None, Dmax=None, pad=1.):
	"""
	"""
	fov_Lx, fov_Ly, boxsize =  get_fov(sample_dimensions, optics_params, detector_params, pdb_file=pdb_file, Dmax=Dmax, pad=pad)
	#
	# deduce number of particle in each direction
	fov_Nx = np.floor(fov_Lx/boxsize)
	fov_Ny = np.floor(fov_Ly/boxsize)
	# 
	x_origin   = -fov_Lx/2. + boxsize/2.
	x_frontier = x_origin + fov_Nx*boxsize
	y_origin   = -fov_Ly/2. + boxsize/2.
	y_frontier = y_origin + fov_Ny*boxsize
	#
	x_range = np.arange(x_origin, x_frontier, boxsize)
	y_range = np.arange(y_origin, y_frontier, boxsize)
	n_particles = np.int(fov_Nx*fov_Ny)
	return x_range, y_range, n_particles
#
def microgaph2particles(micrograph, sample_dimensions, optics_params, detector_params, pdb_file=None, Dmax=30, pad=5.):
	"""
	"""
	fov_Lx, fov_Ly, boxsize =  get_fov(sample_dimensions, optics_params, detector_params, pdb_file=pdb_file, Dmax=Dmax, pad=pad)
	fov_Nx = np.floor(fov_Lx/boxsize)
	fov_Ny = np.floor(fov_Ly/boxsize)
	#
	pixel_size = (fov_Lx/micrograph.shape[1] + fov_Ly/micrograph.shape[0])/2.
	n_boxsize = np.int(boxsize/pixel_size)
	Nx = np.int(fov_Nx*n_boxsize)
	Ny = np.int(fov_Ny*n_boxsize)
	data = micrograph[0:Ny, 0:Nx]
	#print(fov_Lx, fov_Ly, micrograph.)
	particles = slicenstack(data, n_boxsize=n_boxsize)
	return particles
#
def slicenstack(data, n_boxsize=256, n_ovl=0):
    """ slicenstack: reads a 2D numpy array and returns a 3D numpy array
    """
    if n_ovl == 0:
        data_stack = blockshaped(data, n_boxsize, n_boxsize)
    else:
        n_split = math.floor((data.shape[0]-2*n_ovl)/(n_boxsize))
        n_dilat = n_boxsize + 2*n_ovl
        data_stack = np.zeros((n_split*n_split,n_dilat,n_dilat))
        print("Array dimensions: ",data_stack.shape)
        i_stack = 0
        for i in np.arange(n_split):
            for j in np.arange(n_split):
                istart = i*n_boxsize
                istop  = istart + n_dilat
                jstart = j*n_boxsize
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
        n_boxsize = stack.shape[1]
        n = n_split*n_boxsize
        data_unstacked = unblockshaped(stack, n, n)
    else:
        n_boxsize = stack.shape[1] - 2*n_ovl
        n = n_split*n_boxsize
        stacks  = np.arange(0,stack.shape[0])
        rows    = np.arange(n_ovl,n_ovl+n_boxsize)
        columns = np.arange(n_ovl,n_ovl+n_boxsize)
        stack_new = stack[np.ix_(stacks,rows,columns)]
        print(stack_new.shape)
        data_unstacked = unblockshaped(stack_new, n, n)
    print("Array dimensions: ",data_unstacked.shape)
    return data_unstacked	
#
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
#
def get_fov(sample_dimensions, optics_params, detector_params, pdb_file=None, Dmax=None, pad=1.):
	"""
	"""
	# retrieve arguments (lenths in nm)
	detector_Nx           = detector_params[0]
	detector_Ny           = detector_params[1]
	detector_pixel_size   = detector_params[2]*1e3
	magnification         = optics_params[0]
	hole_diameter         = sample_dimensions[0]
	hole_thickness_center = sample_dimensions[1]
	hole_thickness_edge   = sample_dimensions[2]
	# define physical size of field-of-view (in nm)
	detector_Lx = detector_Nx*detector_pixel_size
	detector_Ly = detector_Ny*detector_pixel_size
	fov_Lx = detector_Lx/magnification
	fov_Ly = detector_Ly/magnification
	# define particle boxsize (in nm)
	if Dmax is None:
		if pdb_file is not None:
			Dmax = get_Dmax(pdbfile)
		else:
			Dmax = 100
	boxsize = Dmax + 2*pad
	#
	return fov_Lx, fov_Ly, boxsize
#
def define_grid(sample_dimensions,pdb_file=None,Dmax=None):
    """ define_grid
    """
    if Dmax is None:
        if pdb_file is not None:
            Dmax = get_Dmax(pdbfile)
        else:
            Dmax = 100
    print("Maximum dimension of the molecule: {0} nm".format(Dmax))
    step_size = int(Dmax*3/2)
    xy_origin = int(sample_dimensions[0]/(2*np.sqrt(2)))
    xy_range  = np.arange(-xy_origin, xy_origin+step_size, step_size)
    numpart   = 4*int(((xy_origin+1)/step_size)**2)
    print("Grid of {0} molecules that cover the volume: {1} nm --> {2} nm (step size: {3} nm)".format(numpart,-xy_origin,xy_origin,step_size))
    return xy_range, numpart
#
