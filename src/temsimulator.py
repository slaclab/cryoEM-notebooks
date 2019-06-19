#!/usr/local/bin/python

import os
import sys
import subprocess
import numpy as np
import math
from scipy.stats import special_ortho_group

numpart = 500
pdb = '5t2b'
xandy = 40
rota = 50
dose = 1000

# Checks if a matrix is a valid rotation matrix.
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

for defoc in np.arange(1,3,0.5):

  rotlist = []
  for x in range(0,numpart+1):
	x = special_ortho_group.rvs(3)
	y = rotationMatrixToEulerAngles(x)
	rotlist.append(y)

  for rot in range(0,rota+1,5):
     print defoc, rot
     if not os.path.exists("%s_%s_%s_nr.mrc"%(pdb, rot, defoc)):
	outp = open("%s_%s_%s_nr.txt"%(pdb, rot, defoc), "w")

	outp.write('# File created by TEM-simulator, version 1.3.\n')
	outp.write('%s  6\n'%numpart)
	outp.write('#            x             y             z           phi         theta           psi  \n')

	l = 0
	x = -100
  	for y in range(-400, 500, xandy):
		for x in range(-400, 500, xandy):
				if l == int(numpart):
					break
				outp.write('%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n'%(x, y, 0, rotlist[l][0], rotlist[l][1], rotlist[l][2]))
				l += 1	
	outp.close()
	
	with open('apof.ini', 'r') as file:
		ai = file.read()

	ai = ai.replace("num_particles = 1000", "num_particles = %s"%numpart)
	ai = ai.replace("pdb_file_in = 4v1w.pdb", "pdb_file_in = %s.pdb"%(pdb))
	ai = ai.replace("coord_file_in = 100mol.txt", "coord_file_in = %s_%s_%s_nr.txt"%(pdb, rot, defoc))
	ai = ai.replace("image_file_out = highdose.mrc", "image_file_out = %s_%s_%s_nr.mrc"%(pdb, rot, defoc))
	ai = ai.replace("defocus_nominal = 1", "defocus_nominal = %s"%(defoc))
	ai = ai.replace("dose_per_im = 1000000", "dose_per_im = %s"%(dose))

	with open("outp_%s_%s_%s_nr.ini"%(pdb,rot, defoc), 'w') as file1:
  		file1.write(ai)
	file1.close()

	os.system("~/app/TEM-simulator_1.3/src/TEM-simulator outp_%s_%s_%s_nr.ini"%(pdb, rot, defoc))

     if not os.path.exists("%s_%s_%s_r.mrc"%(pdb, rot, defoc)):
	outp = open("%s_%s_%s_r.txt"%(pdb, rot, defoc), "w")

	outp.write('# File created by TEM-simulator, version 1.3.\n')
	outp.write('%s  6\n'%numpart)
	outp.write('#            x             y             z           phi         theta           psi  \n')

	l2 = 0
	x2 = -100
  	for y2 in range(-400, 500, xandy):
		for x2 in range(-400, 500, xandy):
				if l2 == int(numpart):
					break
				outp.write('%14.4f%14.4f%14.4f%14.4f%14.4f%14.4f\n'%(x2, y2, 0, rotlist[l2][0], float(rotlist[l2][1])+(float(rot)/10), rotlist[l2][2]))
				l2 += 1	
	outp.close()
	
	with open('apof.ini', 'r') as file:
		ai = file.read()

	ai = ai.replace("num_particles = 1000", "num_particles = %s"%numpart)
	ai = ai.replace("pdb_file_in = 4v1w.pdb", "pdb_file_in = %s.pdb"%(pdb))
	ai = ai.replace("coord_file_in = 100mol.txt", "coord_file_in = %s_%s_%s_r.txt"%(pdb, rot, defoc))
	ai = ai.replace("image_file_out = highdose.mrc", "image_file_out = %s_%s_%s_r.mrc"%(pdb, rot, defoc))
	ai = ai.replace("defocus_nominal = 1", "defocus_nominal = %s"%(defoc))
	ai = ai.replace("dose_per_im = 1000000", "dose_per_im = %s"%(dose))

	with open("outp_%s_%s_%s_r.ini"%(pdb,rot, defoc), 'w') as file1:
  		file1.write(ai)
	file1.close()

	os.system("~/app/TEM-simulator_1.3/src/TEM-simulator outp_%s_%s_%s_r.ini"%(pdb, rot, defoc))

for defoc in np.arange(1,3,0.5):
     if not os.path.exists("%s_%s_nr.mrcs"%(pdb, defoc)) or not os.path.exists("%s_%s_r.mrcs"%(pdb, defoc)):
	os.system("newstack %s_%s_%s_nr.mrc %s_%s_%s_nr.mrc %s_%s_%s_nr.mrc %s_%s_%s_nr.mrc %s_%s_%s_nr.mrc %s_%s_%s_nr.mrc %s_%s_%s_nr.mrc %s_%s_%s_nr.mrc %s_%s_%s_nr.mrc %s_%s_%s_nr.mrc %s_%s_%s_nr.mrc %s_%s_nonrot.mrcs"%(pdb, 0, defoc, pdb, 5, defoc, pdb, 10, defoc, pdb, 15, defoc, pdb, 20, defoc, pdb, 25, defoc, pdb, 30, defoc, pdb, 35, defoc, pdb, 40, defoc, pdb, 45, defoc, pdb, 50, defoc, pdb, defoc))
	os.system("newstack %s_%s_%s_r.mrc %s_%s_%s_r.mrc %s_%s_%s_r.mrc %s_%s_%s_r.mrc %s_%s_%s_r.mrc %s_%s_%s_r.mrc %s_%s_%s_r.mrc %s_%s_%s_r.mrc %s_%s_%s_r.mrc %s_%s_%s_r.mrc %s_%s_%s_r.mrc %s_%s_rot.mrcs"%(pdb, 0, defoc, pdb, 5, defoc, pdb, 10, defoc, pdb, 15, defoc, pdb, 20, defoc, pdb, 25, defoc, pdb, 30, defoc, pdb, 35, defoc, pdb, 40, defoc, pdb, 45, defoc, pdb, 50, defoc, pdb, defoc))
