#
import os
import CifFile
from CifFile import ReadCif
#
def star_reader(filename):
    """
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


