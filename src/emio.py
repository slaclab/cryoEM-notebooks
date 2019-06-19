def simio(pdbdir, pdb_keyword, crd_keyword, mrcdir, mrc_keyword):
    """ simio
    """
    pdb_file = pdbdir+pdb_keyword+'.pdb'
    crd_file = mrcdir+crd_keyword+'.txt'
    mrc_file = mrcdir+mrc_keyword+'.mrc'
    log_file = mrcdir+mrc_keyword+'.log'
    inp_file = mrcdir+mrc_keyword+'.inp'
#    
    return pdb_file, mrc_file, crd_file, log_file, inp_file
