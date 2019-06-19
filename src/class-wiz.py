#### Script to investigate convergence behaviour in relion 3D classification
#### Copyright Cornelius Gati 2018 - SLAC Natl Acc Lab - cgati@stanford.edu

import os
import sys
import numpy as np
import matplotlib
#matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import operator
import collections
import pylab as P
import matplotlib.backends.backend_pdf
from operator import itemgetter

################INPUT

print 'Please specify:'
print '--f 		folder path 						(default: current folder)'
print '--root 		root name 						(default: \'run\')'
print '--o		output pdf name 					(default: output.pdf)'
print '--plot		\'bar\' or \'line\'  					(default: bar)'
print '--filt 		\'true\' or \'false\' obtain filtered.star file 	(default: false)'
print '--mic		minimum resolution for Gctf resolution estimate 	(default: none)'

folder = '.'
rootname = 'run'
output = 'output.pdf'
plottype = 'bar'
micfilt = ''
filtstar = 'false'

for si, s in enumerate(sys.argv):
	if s == '--f':
		folder = sys.argv[si+1]
	
	if s == '--root':
		rootname = sys.argv[si+1]
	
	if s == '--o':
		output = sys.argv[si+1]
	
	if s == '--plot':
		plottype = sys.argv[si+1]

	if s == '--filt':
		filtstar = sys.argv[si+1]

	if s == '--mic':
		micfilt = sys.argv[si+1]


#### List all files in folder and sort by name
filesdir = sorted(os.listdir(folder))
unwanted = [];

##Check number of iterations 
iterationlist = []; iterationtemp = []; iterations = []; 
filesdir = sorted(filesdir)

for datafile in filesdir:
	if 'data.star' in datafile and 'sub' not in datafile:
	   if 'it001' in datafile and len(folder) == 0:	#Get rootname by looking at first iteration if not specified
	   	rootname = datafile.split('_it')[0]
	   if 'ct' in datafile.split('_')[-3]:
		if int(datafile.split('_')[-2][2:]) != int(datafile.split('_')[-3][2:]):
			iterationtemp.append(datafile)
	   if 'ct' not in datafile.split('_')[-3]:
			iterationtemp.append(datafile)
	   #if 'ct' in datafile.split('_')[-3]:
		#if int(datafile.split('_')[-2][2:]) != int(datafile.split('_')[-3][2:]):
		#	iterationtemp.append(datafile)
	   iterations.append(int(datafile.split('_')[-2][2:]))
			
### Open PDF for output of figures
pdf = matplotlib.backends.backend_pdf.PdfPages('%s'%output)

## List of all input files used
iterationlist = sorted(iterationtemp, key = lambda x: x.split('_')[-2])
if len(iterations) == 0:
	print ''
	print 'I cannot find any data.star files in the provided folder!'
	sys.exit()

iterations = max(iterations)+1

##Print used input data.star files
print ''
for files in iterationlist:
	print 'Using %s as input'%files

##Check number of particles, number of classes, number of micrographs
part = 0; classes = []; micrographlist = []; checklist = []; checklistcol = []; anglerot = []; rescol = [];

with open('%s/%s_it002_data.star'%(folder, rootname), 'rb') as f:
     	for l in f:
		if 10 < len(l) < 50:  #check header
		   if l.split()[0] == '_rlnClassNumber': #check header for ClassNumber column
			classcolumn = int(l.split()[1][1:])-1		 
		   if l.split()[0] == '_rlnMicrographName': #check header for MicrographName column
			miccolumn = int(l.split()[1][1:])-1
		   if l.split()[0] == '_rlnImageName': #check header for ParticleName column
			particolumn = int(l.split()[1][1:])-1
		   if l.split()[0] == '_rlnAngleRot': #check header for ParticleName column
			anglerot = int(l.split()[1][1:])-1
 		   if l.split()[0] == '_rlnCtfMaxResolution': #check header for CtfMaxResolution column
			rescol = int(l.split()[1][1:])-1

		   if '_rln' in l.split()[0]: 
			checklist.append(int(l.split()[1][1:])-1)
			checklistcol.append(l.split()[0])

	 	if len(l) >= 50:
			part+=1
			classes.append(int(l.split()[classcolumn]))
			micrographlist.append(l.split()[miccolumn]) 

classes = max(classes)

print ''
print 'Plots will be generated for the following columns:', checklistcol


## Initial colorbar
fig = plt.figure()
ax1 = fig.add_axes([0.05, 0.80, 0.9, 0.15])

cmap = plt.get_cmap('jet', int(classes))
norm = matplotlib.colors.Normalize(vmin=1, vmax=int(classes)+1)
ticks = np.arange(1, int(classes)+1)
cb1 = matplotlib.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='horizontal', ticks=ticks)
cb1.set_label('Color for each class')

pdf.savefig()
#plt.show()


###### Go into each iteration_data.star file and read in information, such as particle class assignments etc.
groupnumarray = np.zeros((part, iterations), dtype=np.double) # Class assignments over all iterations
checkarray = np.empty((part, len(checklist)), dtype=np.double) # Stats of final iteration
check = np.empty((part, 2, iterations), dtype=np.double) # Stats of final iteration
changes = []; 

micnumdict = collections.defaultdict(list)

print ''
for datafile in iterationlist:
	miciterdict = collections.defaultdict(list)
	if 'ct' in datafile:
	   	iteration = int(datafile.split('_')[-2][2:])
	if 'ct' not in datafile:
		iteration = int(datafile.split('_')[-2][2:])
   	particle = 0
	changesum = 0
	if int(iteration) > 1:
   		  with open('%s/%s'%(folder, datafile), 'rb') as f:
     			for l in f:
	 			if len(l) >= 50:
					groupnum = l.split()[classcolumn]	# Class number
					micrograph = l.split()[miccolumn]  	## Micrograph name
					particlename = l.split()[particolumn] 	## Particle Name

					groupnumarray[particle, iteration] = groupnum
					miciterdict[micrograph].append(groupnum)
					if int(iteration) == iterations-1: 	#last iteration: all columns of star file to checkarray
					    for i in range(0, len(checklist)-1):
						if 'mrc' not in l.split()[i]:	#No filenames in checklist
							colval = l.split()[i]	
						   	if 'group' in l.split()[i]:
								colval = int(l.split()[i][-2:])
							checkarray[particle, i] = colval
					
					if int(iteration) >= iterations-2:
					    for i in range(0, len(checklist)-1):
							check[particle, 0, iteration] = l.split()[anglerot]
							check[particle, 1, iteration] = l.split()[classcolumn]

					if groupnumarray[particle, iteration] == groupnumarray[particle, int(iteration)-1]:
						change = 'nan'; 
						change1 = 0;
					if groupnumarray[particle, iteration] != groupnumarray[particle, int(iteration)-1]:
						change = groupnum
						change1 = 1;
						changesum += 1;
					
					particle += 1
	changes.append(changesum)

	print "Iteration %s: %s particles changed class assignments"%(iteration, changesum)
	
	## After each iteration create Histogram of class assignments for each micrograph
	miciterdict = collections.OrderedDict(sorted(miciterdict.items(), key=operator.itemgetter(0)))
	for key, value in miciterdict.iteritems():
		michisto = []
		for numb in value:
			michisto.append(float(numb))
		classes = int(classes)
		michisto = np.histogram(michisto, bins=classes, range=(1, classes+1))
		micnumdict[key].append(michisto[0].tolist())

######## Plot rotational and translational accuracy over each iteration

rotationcol = 2; translationcol = 3;
rotation = np.zeros((int(classes)+1, iterations), dtype=np.double); 
translation = np.zeros((int(classes)+1, iterations), dtype=np.double);

for datafile in iterationlist:
	iteration = int(datafile.split('_')[-2][2:])
	with open('%s/%s_model.star'%(folder, datafile[:-10]), 'rb') as f:
		for l in f:
			if '_rlnAccuracyRotations' in l:
				rotationcol = int(l.split('#')[-1])-1
			if '_rlnAccuracyTranslations' in l:
				translationcol = int(l.split('#')[-1])-1
			if 'class' in l and 'mrc' in l and 'mrcs' not in l:
				classnum = int(l.split('.mrc')[0][-3:])
				rotation[classnum, iteration] = l.split()[rotationcol]
				translation[classnum, iteration] = l.split()[translationcol]

rotation = np.array(rotation[1:])
translation = np.array(translation[1:])


if len(set(rotation[0])) == 1:
	print 'You did not perform image alignment during classification - skipping these two plots!'

if len(set(rotation[0])) > 1:

#Rotational

	cmap = plt.get_cmap('jet', int(classes))
	plt.figure(num=None, dpi=80, facecolor='white')
	plt.title('RotationalAccuracy', fontsize=16, fontweight='bold')
	plt.xlabel('Iteration #', fontsize=13)
	plt.ylabel('RotationalAccuracy', fontsize=13)
	plt.grid()
	colors = np.arange(0, int(classes))
	d = 0;
	for c, r in zip(colors, rotation):
		d = c+1
		plt.plot(r[2:], linewidth=3, color=cmap(c), label='Class %s'%d)
	plt.legend(loc='best')
	pdf.savefig()
	#plt.show()

#Translational
	cmap = plt.get_cmap('jet', int(classes))
	plt.figure(num=None, dpi=80, facecolor='white')
	plt.title('TranslationalAccuracy', fontsize=16, fontweight='bold')
	plt.xlabel('Iteration #', fontsize=13)
	plt.ylabel('TranslationalAccuracy', fontsize=13)
	plt.grid()
	colors = np.arange(0, int(classes))
	d=0
	for c, t in zip(colors, translation):
		d = c+1;
		plt.plot(t[2:], linewidth=3, color=cmap(c), label='Class %s'%d)
	plt.legend(loc='best')
	pdf.savefig()
	#plt.show()


###########################################################################

### Sort group assignment array column by column
sortindices = np.lexsort(groupnumarray[:,1:].T)
groupnumarraysorted = groupnumarray[sortindices]

### Heat map of group sizes
H = groupnumarraysorted[:,2:]
cmap = plt.get_cmap('jet', int(classes))
norm = matplotlib.colors.Normalize(vmin=1, vmax=int(classes)+1)
plt.figure(num=None, dpi=120, facecolor='white')
plt.title('Class assignments of each particle', fontsize=16, fontweight='bold')
plt.xlabel('Iteration #', fontsize=13)
plt.ylabel('Particle #', fontsize=13)
mat = plt.imshow(H,aspect='auto', interpolation="nearest", cmap=cmap, norm=norm)
ticks = np.arange(2, iterations+1)
plt.gca().set_xticklabels(ticks)
cb1 = plt.colorbar(mat, ticks=np.arange(1, classes+1))
cb1.set_label('Class #')
pdf.savefig()
#plt.show()

### Plot heat map of the last 5 iterations (close-up)
H = groupnumarraysorted[:,-5:]
cmap = plt.get_cmap('jet', int(classes))
norm = matplotlib.colors.Normalize(vmin=1, vmax=int(classes)+1)
plt.figure(num=None, dpi=120, facecolor='white')
plt.title('Class assignments of each particle - last 5 iterations', fontsize=16, fontweight='bold')
plt.xlabel('Iteration #', fontsize=13)
plt.ylabel('Particle #', fontsize=13)
mat = plt.imshow(H, aspect='auto', interpolation="nearest", cmap=cmap, norm=norm)
ticks = np.arange(iterations-5, iterations+1)
plt.gca().set_xticklabels(ticks)
cb1 = plt.colorbar(mat, ticks=np.arange(1, classes+1))
cb1.set_label('Class #')
pdf.savefig()
#plt.show()

#########################################################################################################################################################
### Sort group assignment array column by column
check = check[sortindices]
#### Jumper analysis
checkdict = collections.defaultdict(list)
checktest = [];

for c in check:
	checkdict[c[:][1][-1]].append(c[:][1][-2])
for key, value in checkdict.iteritems():
	hist, bins = np.histogram(value, bins=np.arange(1, int(classes)+2), normed=True)
	checktest.append(hist)

fig = plt.figure(num=None, dpi=80, facecolor='white')
ax = fig.add_subplot(111)
plt.title('Class assignment of each particle - last iteration', fontsize=16, fontweight='bold')
plt.xlabel('Class assignment iteration %s'%(int(iterations)-1), fontsize=13)
plt.ylabel('Class assignment iteration %s'%(int(iterations)), fontsize=13)
plt.grid()
ticks = np.arange(0, int(classes)+1)
labels = np.arange(1, int(classes)+2)
plt.xticks(ticks, labels)
plt.yticks(ticks, labels)
plt.imshow(checktest, aspect='auto', interpolation="nearest", origin='lower')	#FIXME values in box
cb2 = plt.colorbar(ticks=np.arange(0, 1, 0.1))
cb2.set_label('Fraction of particles went into group #')
plt.figtext(0, 0, 'Particle class assignment changes in the last iteration')
pdf.savefig()
#plt.show()

######################################################################################################################

###########################################################################
#### Class assignments per micrograph of the last iteration

### Convert dictionary to 2D array
micval = []; micfirst = [];
micticks = []
for mics in micnumdict:
	max_index, max_value = max(enumerate(micnumdict[mics][-1]), key=operator.itemgetter(1))
	micval.append(micnumdict[mics][-1])
	micticks.append(mics)

### Sort last iter by most prominent group
micval = np.array(micval)

### Plot heat map last iteration
cmap = plt.get_cmap('jet', np.max(micval)-np.min(micval)+1)
plt.figure(num=None, dpi=80, facecolor='white')
plt.title('Class assignments of each micrograph - last iteration', fontsize=16, fontweight='bold')
plt.xlabel('Class #', fontsize=13)
plt.ylabel('Micrograph #', fontsize=13)
ticks = np.arange(0, int(classes)+1)
labels = np.arange(1, int(classes)+2)
plt.xticks(ticks, labels)
plt.grid()
plt.imshow(micval, aspect='auto', interpolation="nearest", cmap=cmap)
cb3 = plt.colorbar()
cb3.set_label('Total number of particles in class #')
plt.figtext(0, 0, 'Micrographs contributing to certain classes (e.g. important when merging datasets)')
pdf.savefig()
#plt.show()

###### Find out how often particles are jumping
scorelist = [];
for g in groupnumarraysorted[:,2:]:

	count = 1;
	stayed = 0;
	g = g[::-1]  ## reverse order 
	for i,ig in enumerate(g):
		if i > 0 and count > 0:

			if len(set(g)) == 1:
				stayed = len(g);
				count = 0;
			if int(ig) == int(g[i-1]):
				count += 1;
			if int(ig) != int(g[i-1]):
				stayed = count;
				count = 0;
	   
	g = g.tolist()
	counter = [[x,g.count(x)] for x in set(g)]
	#score3 = (iterations-2)/float(len(counter))
	score2 = (float(len(counter)))
	score1 = (stayed/score2)/(iterations-2)
	scorelist.append(score1)

	#print g[::-1], counter, stayed, 'Score1', score1 #, 'Score2', score2#, 'Score3', score3

plt.figure(num=None, dpi=80, facecolor='white')
plt.title('Particle jump score', fontsize=16, fontweight='bold')
plt.xlabel('Score = (# class assignments/# of iterations)', fontsize=13)
plt.ylabel('# of particles with score normalized', fontsize=13)
plt.grid()
#histbins = sorted(set(scorelist))
histbins = np.arange(0, 0.5, 0.05)
plt.hist(scorelist, bins=histbins, normed=True)
mean1 = np.mean(scorelist)
variance1 = np.var(scorelist)
sigma1 = np.sqrt(variance1)
#x1 = np.linspace(min(scorelist), max(scorelist), 100)
x1 = np.linspace(0, 0.5, 100)
plt.figtext(0, 0, 'Gaussian sigma: %s, variance: %s, mean: %s'%(sigma1, variance1, mean1))
plt.plot(x1,mlab.normpdf(x1, mean1, sigma1))
pdf.savefig()

########################
### Throw out particles that change assignments more often than (mu + 2*sigma)
cutoff = (mean1+(2*sigma1))

for i, score in enumerate(scorelist):
	if score > cutoff:
		unwanted.append(i)
########################

### Number of assignment changes per iteration
plt.figure(num=None, dpi=80, facecolor='white')
plt.title('Total assignment changes per iteration', fontsize=16, fontweight='bold')
plt.xlabel('Iteration #', fontsize=13)
plt.ylabel('# of changed assignments', fontsize=13)
plt.grid()
plt.plot(changes[2:])
pdf.savefig()
#plt.show()

#########################################################################################################################################################
checkarraysorted = checkarray[sortindices]
### Plot histogram of each column in data.star sorted by the class assignments of the last iteration
histocolarray = collections.defaultdict(list)
for column in checkarraysorted.transpose():
	lastitgrouparray = collections.defaultdict(list)
	for ci, col in enumerate(column):
		lastitgrouparray[groupnumarraysorted[:,-1][ci]-1].append(col)
	for key, value in lastitgrouparray.iteritems():
		histocolarray[key].append(value)

fincolarray = collections.defaultdict(list)
for key, value in histocolarray.iteritems():
	for vi, v in enumerate(value):
		fincolarray[vi].append(v)

for key2, value2 in fincolarray.iteritems():
	rangeval = []; temp = []; temp2 = [];
	for vi2, v2 in enumerate(value2):
		for v3 in v2:
			rangeval.append(v3)
	if len(set(rangeval)) > 1:
	    for vi2, v2 in enumerate(value2):
		temp.append(v2)
	
	if len(temp) > 0:
		plt.figure(num=None, dpi=120, facecolor='white')
		#####
		
		if plottype == 'line':
		   for ti, t in enumerate(temp):

			n, bins = np.histogram(t, bins=15, range=(min(rangeval),max(rangeval)), normed=True)
			ti += 1;
			plt.plot(bins[:-1], n, linewidth=3, color=cmap(ti), label='Class %s'%ti)
			plt.legend(loc='best')
		#####
		
		if plottype == 'bar':
			n, bins, patches = plt.hist(temp, bins=15, range=(min(rangeval),max(rangeval)), normed=1, histtype='bar')
		
			cmap = plt.get_cmap('jet', int(classes))
			colors = np.arange(0, int(classes))
			for c, p in zip(colors, patches):
			    d = c+1;
			    plt.setp(p, 'facecolor', cmap(c))
			#plt.colorbar(ticks=np.arange(1, int(classes)+1))

		plt.title('Histogram Column %s'%checklistcol[key2], fontsize=16, fontweight='bold')
		plt.xlabel('%s'%checklistcol[key2], fontsize=13)
		plt.ylabel('# particles per bin', fontsize=13)
		plt.grid()

		pdf.savefig()
		#plt.show()

print 'Saved all plots in %s'%output
pdf.close()

################################################################### DELETE UNWANTED PARTICLES FROM INITIAL STAR FILE

particcount = 0;
initstarfile = '%s/%s_it001_data.star'%(folder,rootname)
a1 = open('%s_filtered.star'%(rootname), 'w')
if filtstar != 'false' or micfilt != '':
 with open(initstarfile, 'rb') as g:
   for m in g:
	if len(m) < 50:
		a1.write(m)
	if len(m) > 50:
		if micfilt != '':
			if float(m.split()[rescol]) <= float(micfilt):
				a1.write(m)
				print 'micfilt'
		if filtstar != 'false':
			if particcount not in unwanted:
				a1.write(m)
				print 'nomicfilt'
		if filtstar != 'false' and micfilt != '':
			if float(m.split()[rescol]) <= float(micfilt) and particcount not in unwanted:
				a1.write(m)
				print 'both'

		particcount += 1;
 print 'Saved %s_filtered.star file ommitting %s out of %s particles that changed classes too often'%(rootname, len(unwanted), particcount)
a1.close()

