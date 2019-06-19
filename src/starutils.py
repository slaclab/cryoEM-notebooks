"""
STARUTILS - F.Poitevin, Stanford, 2019
a collection of python tools to manipulate STAR files
"""
import re
import numpy as np
#
def load_star(filename,usecols=None):
    """load_star: we read a star file and only keep the usecols columns in particle lines.
    """
    lines = (line for line in open(filename, 'rb')  if re.search("@", line.decode('utf-8')) )
    data = np.genfromtxt(lines,usecols = usecols,dtype='U')
    return data
