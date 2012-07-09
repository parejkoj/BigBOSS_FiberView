#!/usr/bin/env python
"""
Read in Will Emmett's .txt files, containing position measurements of the 
fibers positions in the invar plate, and compare the measured positions.
"""
# Initial version: 2012.06.20 John Parejko
import csv
import numpy as np
import numpy.core.records as rec
import numpy.lib.recfunctions as recfunc
import pyfits

import csv_parse
import fvmapping
from kdtree_test.naive_pairs import xymatch

import matplotlib
from matplotlib import rc
rc('text', usetex=True)
import matplotlib.pyplot as plt
plt.rcParams.update({'backend':'pdf',
                     'text.usetex':True,
                     'text.fontsize':18,
                     'axes.labelsize': 18,
                     'legend.fontsize': 14,
                     'xtick.labelsize': 18,
                     'ytick.labelsize': 18,
                     'axes.linewidth':2})
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
import matplotlib.gridspec as gridspec
plt.ion() # interactive mode on

# these names and types should be the same on every file.
names = ['step','feature','x','y','z','size']
types = ['int','string','float','float','float','float']

def read_originals():
    """Return the originally defined fiber positions, sorted by x and y."""
    head,points1 = csv_parse.read('../Baltay-fibers_random.csv',delimiter=' ')
    head,points2 = csv_parse.read('../Baltay-fibers_residual.csv',delimiter=' ')
    points2 = recfunc.drop_fields(points2,('r','theta'))
    points2['Number'] += 10000 # to distinguish them from the "randoms"
    return np.hstack((points1,points2))
#...

def read_text(filename):
    """
    Read in .TXT filename and return a numpy record array.
    """
    parser = csv_parse.csv_parse(names,types)
    infile = csv.reader(open(filename,'r'),delimiter=' ',skipinitialspace=True)
    # skip three header lines, use the next one
    for i in range(4):
        line = infile.next()
    data = []
    try:
        while line:
            data.append(parser.convert2types(line))
            line = infile.next()
    except StopIteration:
        # we reached the end of the file
        pass
    return rec.fromrecords(data,names=names,formats=parser._formats,shape=len(data))
#...

def read_fits(filename):
    """Read (FITS) filename and return an array of objects most likely to be
    detected fibers."""
    # NOTE; the images are 8176x6132
    data = pyfits.open(filename)[1].data
    # this limits it to sources that aren't one or two high pixels
    test = (data.A_IMAGE * data.B_IMAGE) > 2.
    test &= data.FLUX_ISO > 2e4
    return data[test]
#...

def coalign_fits(x0,y0,x1,y1):
    """
    Return x1,y1 reordered to be in the same order as x0,y0.
    """    
    match = xymatch(np.array(zip(x0,y0)),np.array(zip(x1,y1)),40,maxmatch=1)
    # Have to keep everything aligned to the original x0,y0 positions.
    match.sort(order='idx1')
    return x1[match['idx2']],y1[match['idx2']]
#...

def make_arrays(data,xname='x',yname='y'):
    """Returns arrays of x and y positions."""
    # they all have the same length
    names = sorted(data.keys())
    x = np.zeros((len(data),len(data[names[0]])))
    y = np.zeros((len(data),len(data[names[0]])))
    for i,n in enumerate(names):
        if i != 0 and 'IMAGE' in xname:
            x[i],y[i] = coalign_fits(x[0],y[0],data[n][xname],data[n][yname])
        else:
            x[i] = data[n][xname]
            y[i] = data[n][yname]

    # NOTE: the images are flipped around the y axis
    if 'IMAGE' in xname:
        y = 8176-y
    return x,y,names
#...

def plot_sep_vector(xy0,xy1,name):
    """Plot vectors from xy1 (drilled) to the xy0 (specified) positions,
    as given by T(xy1)."""
    transform = fvmapping.Transform(xy1,xy0)
    transform()
    points = transform.transform_all()
    match = transform.check()
    diff0 = points[match['idx1']] - xy0[match['idx2']]
    diff1 = xy1[match['idx1']] - points[match['idx1']]
    diff01 = xy1[match['idx1']] - xy0[match['idx2']]
    r = np.sqrt((xy0[match['idx2']]**2).sum(axis=1))
    rdiff0 = np.sqrt((diff0**2).sum(axis=1))
    rdiff1 = np.sqrt((diff1**2).sum(axis=1))
    rdiff01 = np.sqrt((diff01**2).sum(axis=1))

    extent = (-100,600,-150,300)
    xticks = np.arange(-100,600,100)
    fig = plt.figure(figsize=(16,8))
    gs = gridspec.GridSpec(1,2,wspace=0)
    ax_vec = fig.add_subplot(gs[0])
    plt.xticks(xticks)
    ax_resid = fig.add_subplot(gs[1],sharey=ax_vec)
    ax_vec.scatter(xy0[match['idx2'],0],xy0[match['idx2'],1],s=4,c='black',edgecolors='none')
    # NOTE: quiver plots arrows from x,y of length u,v, not from x,y to u,v!
    scale = 2e-3
    Q1 = ax_vec.quiver(xy0[match['idx2'],0],xy0[match['idx2'],1],diff01[:,0],diff01[:,1],color='green',width=0.005,label='coords',scale=scale,scale_units='xy')
    Q2 = ax_vec.quiver(xy1[match['idx1'],0],xy1[match['idx1'],1],diff1[:,0],diff1[:,1],color='blue',width=0.005,label='transform',scale=scale,scale_units='xy')
    # the physical length is U/scale: 0.1/2e-3 = 50 microns
    ax_vec.quiverkey(Q1,0.1,0.97,0.1,r'$50\,\mu m$',labelpos='E',fontproperties={'size':14})
    ax_vec.quiverkey(Q2,0.1,0.92,0.1,r'$50\,\mu m$',labelpos='E',fontproperties={'size':14})
    ax_vec.set_title('coordinate difference and transformed difference')
    ax_vec.set_xlabel('x (mm)')
    ax_vec.set_ylabel('y (mm)')
    ax_vec.axis(extent)
    Q3 = ax_resid.quiver(xy0[match['idx2'],0],xy0[match['idx2'],1],diff0[:,0],diff0[:,1],color='red',width=0.005,label='residual',scale=scale,scale_units='xy')
    ax_resid.quiverkey(Q3,0.1,0.97,0.1,r'$50\,\mu m$',labelpos='E',fontproperties={'size':14})
    ax_resid.set_title('post-transform residuals')
    ax_resid.set_xlabel('x (mm)')
    ax_resid.axis(extent)
    # TBD: don't really need this, right?
    #fig.suptitle(r'specified $\rightarrow$ '+name,fontsize=18)
    fig.subplots_adjust(hspace=0)
    plt.setp([a.get_yticklabels() for a in fig.axes[1:]], visible=False)
    filename = '../plots/separation_vector-'+name+'.png'
    plt.savefig(filename,bbox_inches='tight',dpi=150)
    print 'Generated:',filename

    # now some histograms
    fig = plt.figure()
    ax_hist = fig.add_subplot(111)
    x_range = (0,0.35)
    ax_hist.hist(rdiff1,bins=50,lw=2,label='coords',histtype='step',range=x_range)
    ax_hist.hist(rdiff01,bins=50,lw=2,label='transform',histtype='step',range=x_range)
    ax_hist.hist(rdiff0,bins=50,lw=2,label='residual',histtype='step',range=x_range)
    plt.legend(loc='upper right')
    plt.xlabel('r (mm)')
    plt.ylabel('number')
    filename = '../plots/separation_hist-'+name+'.png'
    plt.savefig(filename,bbox_inches='tight')
    print 'Generated:',filename

    # now some line plots
    fig = plt.figure()
    ax_diff = fig.add_subplot(111)

    from scipy.optimize import leastsq
    f = lambda v,x: v[0]*x + v[1]
    e = lambda v,x,y: f(v,x) - y
    v0 = [1,0]
    xfit,success = leastsq(e,v0,args=(xy0[:,0],diff0[:,0]))
    yfit,success = leastsq(e,v0,args=(xy0[:,1],diff0[:,1]))
    rfit,success = leastsq(e,v0,args=(r,rdiff0))

    def plot_fit(x,y,fit,label,color):
        xx = np.arange(x.min(),x.max(),.1)
        ax_diff.scatter(x,y,c=color,s=4,edgecolor='none')
        ax_diff.plot(xx,f(fit,xx),c=color,lw=2,label=label+': %4.2e'%fit[0])
    #...
    plot_fit(xy0[match['idx2'],0],diff0[:,0],xfit,'x','red')
    plot_fit(xy0[match['idx2'],1],diff0[:,1],yfit,'y','green')
    plot_fit(r,rdiff0,rfit,'r','blue')
    plt.xlabel('coordinate (mm)')
    plt.ylabel('offset from specified position (mm)')
    plt.title('residual offsets after transforming: '+name)
    plt.legend(loc='upper left')
    filename = '../plots/separation_offsets-'+name+'.png'
    plt.savefig(filename,bbox_inches='tight')
    print 'Generated:',filename
    #import pdb
    #pdb.set_trace()
#...

def hist_sigma(data,label,N,bins=50,range=range):
    plt.hist(data.std(axis=0),lw=2,bins=bins,histtype='step',label=label,range=range)
#...

def plot_sep_hist(x,y,isFits=False):
    """
    Plot the variation in the measurements.
    Set isfits to force the points to be pre-matched before calculations,
    and to adjust plot labels and filename.
    """
    ind = np.arange(64)+1
    zeros = np.zeros(64)
    N = len(x[:,0])
    r = np.sqrt(x**2+y**2)
    plt.figure()
    if isFits:
        bins = 20
    else:
        bins = 50
    range = (min((x.std(axis=0).min(),y.std(axis=0).min(),r.std(axis=0).min())),
             max((x.std(axis=0).max(),y.std(axis=0).max(),r.std(axis=0).max())))
    hist_sigma(x,'x',N,bins,range)
    hist_sigma(y,'y',N,bins,range)
    hist_sigma(r,r'$r=(x^2+y^2)$',N,bins,range)
    if isFits:
        units = '\mathrm{pixels}'
    else:
        units = '\mathrm{mm}'
    plt.xlim(range)
    plt.xlabel(r'$\sigma_{\mathrm{'+str(N)+'\;measurements}} ('+units+')$')
    plt.ylabel('number')
    plt.title('x,y,r standard deviation')
    plt.legend()
    suffix = str(N)
    if isFits:
        suffix += 'fits_'
    else:
        suffix += 'plate_'
    filename = '../plots/xyr-stddev_'+suffix+'.png'
    plt.savefig(filename)
    print 'Generated:',filename
#...

def plot_image_values(data):
    """Make some histograms of various values from the sextracted images."""
    Nimages = len(data)
    Npoints = len(data[data.keys()[0]])
    flux = np.zeros(Nimages*Npoints,'f8')
    for i,n in enumerate(data):
        start = i*Npoints
        end = (i+1)*(Npoints)
        flux[start:end] = data[n]['FLUX_ISO']

    plt.figure()
    plt.hist(flux,histtype='step',bins=20,lw=2)
    plt.xlabel('flux (total counts): '+str(Nimages)+' measurements')
    plt.ylabel('number')
    plt.ylim(0,100)
    filename = '../plots/flux_hist.png'
    plt.savefig(filename)
    print 'Generated:',filename
#...

def main(argv=None):
    if argv is None: argv = sys.argv[1:]
    from optparse import OptionParser, OptionGroup

    usage = '%prog [OPTIONS] FILES'
    usage += '\nRead in xyz maps of the fiber plane and determine how far apart the measurements are.'
    usage += '\nIf the files end in .txt, assume they are measurements from the plate measuring machine.'
    usage += '\nIf the files end in .fits, assume they are sextractor fits to camera images.'
    parser = OptionParser(usage)
    #parser.add_option('-s','--save',dest='save',action='store_true',default=False,
    #                  help='Save the sextractor output .fits files (%default).')
    (opts,args) = parser.parse_args(argv)

    if len(args) == 0:
        parser.error('Need to pass some measured position files!')

    orig = read_originals()
    xy_orig = np.array(zip(orig['x'],orig['y']))

    data = {}
    if 'txt' in args[0] or 'TXT' in args[0]:
        for f in args:
            name = f.split('pass')[-1].split('good')[0]
            data[name] = read_text(f)
            xname = 'x'
            yname = 'y'
            isFits = False
    elif 'fits' in args[0]:
        for f in args:
            name = f.split('_sex')[0][-2:]
            data[name] = read_fits(f)
            xname = 'X_IMAGE_DBL'
            yname = 'Y_IMAGE_DBL'
            isFits = True
    else:
        parser.error("I don't know how to process files of type: "+os.path.splitext(args[0])[1])

    x,y,names = make_arrays(data,xname=xname,yname=yname)
    plot_sep_hist(x,y,isFits=isFits)
    if isFits:
        plot_image_values(data)

    for i,n in enumerate(names):
        plot_sep_vector(xy_orig,np.array(zip(x[i],y[i])),n)
    import pdb
    pdb.set_trace()
#...

if __name__ == "__main__":
    import sys
    sys.exit(main())
