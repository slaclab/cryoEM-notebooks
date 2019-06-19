import numpy as np
import mdtraj as md
from scipy.spatial.distance import pdist
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
