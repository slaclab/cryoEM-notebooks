#
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator
#
def view_image(img):
	"""
	"""
	fig = plt.figure(figsize=(6,6))
	plt.imshow(img, cmap='Greys')
	plt.colorbar()
#
def view_particles(data, slicing=(1,1,1), figsize=1, ncol=5 ):
	"""
	"""
	view = data[::slicing[0],::slicing[1],::slicing[2]]
	figsize = int(figsize*ncol)
	nrow = np.ceil(view.shape[0]/ncol)
	fig = plt.figure( figsize=(ncol*figsize,nrow*figsize) )
	for i in np.arange( view.shape[0] ):
		fig.add_subplot( nrow, ncol, i+1 )
		plt.imshow(view[i], cmap='Greys')
	plt.tight_layout()
	plt.show()

############
# < BIPLOTS

def biplots(prj,prj2=None,n=1,plottype='hexbin',nbins=10,figsize=-1,c=None,show_histo=False,figname=''):
    """ biplots : plot populations in component space

    Description
    -----------
    For n components, shows biplots between all pairs of prj in upper triangle.
    If prj2 is provided, shows biplots between all pairs of prj2 in lower one.
    The possibility to color based on input assignment is offered.

    """
    if c is not None:
        plottype='scatter'
    if(plottype=='scatter'):
        cmap='rainbow'
    else:
        cmap='plasma'
    if(figsize < 0 ):
        if(n == 1):
            figsize=1
        else:
            figsize=4
    figsize=figsize*6
    labels = get_labels(n)
    nrow=n
    ncol=n
    gs = gridspec.GridSpec(nrow, ncol, hspace=0, wspace=0)
    minorticklocator = MultipleLocator(0.1)
    majorticklocator = MultipleLocator(0.2)
    fig = plt.figure(figsize=(figsize,figsize), dpi= 160, facecolor='w', edgecolor='k')
    nbins_coarse = int(nbins/1)
    for i in np.arange(0,nrow,1):
        for j in np.arange(0,ncol,1):
            if(i<j):
                ax = biplot_axis(n,i,j,gs,fig,labels,Ax=prj[:,j],Ay=prj[:,i],c=c,cmap=cmap,nbins=nbins,majortick=0.2,minortick=0.1,linewidth=2.5,plottype=plottype)
            elif(i==j):
                if(show_histo):
                    ax = fig.add_subplot(gs[i,j])
                    plt.grid()
                    plt.tick_params(axis='both', which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False)
                    Ax = prj[:,i]
                    plt.hist(Ax,bins=nbins_coarse)
                    if prj2 is not None:
                        Ay = prj2[:,i]
                        plt.hist(Ay,bins=nbins_coarse,rwidth=0.4)
                else:
                    if(i==0 and c is not None):
                        ax = fig.add_subplot(gs[i,j])
                        ax.set_xlabel(labels[j],fontsize='xx-large')
                        ax.set_ylabel(labels[i],fontsize='xx-large')
                        ax.xaxis.set_label_position('top')
                        for corner in ('top','bottom','right','left'):
                            ax.spines[corner].set_linewidth(0)
                        ax.tick_params(axis='both', which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False,labelsize='large')
                        colorbar = range(1,np.max(c)+1,1)
                        ax.scatter(colorbar, colorbar, c=colorbar, vmin=1, vmax=np.max(c), cmap=cmap)
                        for idx, txt in enumerate(colorbar):
                            ax.annotate(txt,(colorbar[idx],colorbar[idx]))
            else:
                if prj2 is not None:
                    ax = biplot_axis(n,i,j,gs,fig,labels,Ax=prj2[:,j],Ay=prj2[:,i],c=c,cmap=cmap,nbins=nbins,majortick=0.2,minortick=0.1,linewidth=2.5,plottype=plottype)
    plt.tight_layout()
    plt.show()
    if(figname):
        fig.savefig(figname+'_biplot.png')

def biplot_axis(n,i,j,gs,fig,labels,Ax=np.zeros(1),Ay=np.zeros(1),c=None,cmap=None,nbins=1,majortick=0.2,minortick=0.1,linewidth=2.5,plottype='scatter'):
    minorticklocator = MultipleLocator(minortick)
    majorticklocator = MultipleLocator(majortick)
    ax = fig.add_subplot(gs[i,j])
    ax.xaxis.set_minor_locator(minorticklocator)
    ax.xaxis.set_major_locator(majorticklocator)
    ax.yaxis.set_minor_locator(minorticklocator)
    ax.yaxis.set_major_locator(majorticklocator)
    for corner in ('top','bottom','right','left'):
        ax.spines[corner].set_linewidth(linewidth)
    ax.minorticks_on()
    ax.grid(which='major',axis='both',linestyle='-',linewidth=0.5)
    ax.grid(which='minor',axis='both',linestyle='--',linewidth=0.1)
    ax.axvline(0, linestyle=':', linewidth=2, color='k')
    ax.axhline(0, linestyle=':', linewidth=2, color='k')
    if(i == 0 and j != n-1):
        ax.set_xlabel(labels[j],fontsize='xx-large')
        ax.xaxis.set_label_position('top')
        if(j==1):
            ax.tick_params(axis='both', which='both',bottom=False,top=False,left=True,right=False,labelbottom=False,labeltop=False,labelleft=True,labelright=False,labelsize='large')
        else:
            ax.tick_params(axis='both', which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False,labelsize='large')
    elif(i == 0 and j == n-1):
        ax.set_xlabel(labels[j],fontsize='xx-large')
        ax.xaxis.set_label_position('top')
        ax.tick_params(axis='both', which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False,labelsize='large')
    elif(i != 0 and j == n-1):
        ax.tick_params(axis='both', which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False,labelsize='large')
    elif(i < j and i != 0 and j != n-1):
        ax.tick_params(axis='both', which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False,labelsize='large')
    if(j == 0 and i != n-1 ):
        ax.set_ylabel(labels[i],fontsize='xx-large')
        if(i==1):
            ax.xaxis.set_label_position('top')
            ax.tick_params(axis='both', which='both',bottom=False,top=True,left=False,right=False,labelbottom=False,labeltop=True,labelleft=False,labelright=False,labelsize='large')
        else:
            ax.tick_params(axis='both', which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False,labelsize='large')
    elif(j == 0 and i == n - 1):
        ax.set_ylabel(labels[i],fontsize='xx-large')
        ax.tick_params(axis='both', which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False,labelsize='large')
    elif(j != 0 and i == n-1 ):
        ax.tick_params(axis='both', which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False,labelsize='large')
    elif(j < i and j !=0 and i != n-1):
        ax.tick_params(axis='both', which='both',bottom=False,top=False,left=False,right=False,labelbottom=False,labeltop=False,labelleft=False,labelright=False,labelsize='large')
    if(plottype == 'scatter'):
        ax.scatter(Ax, Ay, c=c, cmap=cmap)
    else:
        ax.hexbin(Ax, Ay, gridsize=nbins, cmap=cmap, mincnt=1)
    return ax

def get_labels(n):
    labels = []
    for i in np.arange(0,n,1):
        labels.append('# '+str(i+1))
    return labels

