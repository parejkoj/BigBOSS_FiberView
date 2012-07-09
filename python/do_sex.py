#!/usr/bin/env python
"""
Run Sextractor on simulated files and compare with known positions.
"""

import subprocess
import os
import os.path
import re
import numpy as np
import numpy.lib.recfunctions as rec
import pylab
import pyfits

import csv_parse

import fvmapping

def run_sex(filename):
    """
    Run sextractor and read-in and return the output.
    Note that sextractor sees the points in pixel xy coordinates,
    and the axes are flipped.
    """
    outfile = '../data/sex/'+os.path.splitext(os.path.basename(filename))[0]+'_sex.fits'
    logfile = '../temp/'+os.path.basename(outfile)
    cmd = ' '.join(('sex',filename,'-CATALOG_NAME',outfile))
    stdout=open(logfile+'.log','w')
    stderr=open(logfile+'.err','w')
    subprocess.call(cmd,shell=True,stdout=stdout,stderr=stderr)
    return np.array(pyfits.open(outfile)[1].data),outfile
#...

def make_regions(measured,sexfilename):
    """Make a ds9 region file from the sextractor output."""
    try:
        x = measured['X_IMAGE_DBL']
        y = measured['Y_IMAGE_DBL']
    except KeyError:
        x = None
        y = None
    try:
        a = measured['A_IMAGE']
        b = measured['B_IMAGE']
        t = measured['THETA_IMAGE']
    except KeyError:
        a = None
        b = None
        t = None
    try:
        xpeak = measured['XPEAK_IMAGE']
        ypeak = measured['YPEAK_IMAGE']
    except KeyError:
        xpeak = None
        ypeak = None
    try:
        flux = measured['FLUX_ISO']
    except KeyError:
        flux = None

    if a is not None:
        outfile = open(sexfilename.replace('.fits','-ellipses.reg'),'w')
        # write some default parameters
        outfile.write('global color=cyan font="helvetica 10 normal" edit=1 move=0 delete=1 include=1\n')
        formatter = lambda x,y,a,b,t,flux: 'ellipse('+','.join((str(x),str(y),str(a),str(b),str(t)))+') #text="%.2e"'%flux
        outfile.write('\n'.join(map(formatter,x,y,a,b,t,flux)))
        outfile.write('\n')
        formatter = lambda x,y: 'point('+','.join((str(x),str(y)))+') #point=cross width=2'
        outfile.write('\n'.join(map(formatter,x,y)))
        outfile.close()
    if xpeak is not None:
        outfile = open(sexfilename.replace('.fits','-peaks.reg'),'w')
        # write some default parameters
        outfile.write('global color=red font="helvetica 12 normal" edit=1 move=0 delete=1 include=1 width=2\n')
        formatter = lambda xpeak,ypeak: 'point('+','.join((str(xpeak),str(ypeak)))+') #point=x'
        outfile.write('\n'.join(map(formatter,xpeak,ypeak)))
        outfile.write('\n')
        outfile.close()
#...

def read_positions():
    head,points1 = csv_parse.read('../Baltay-fibers_random.csv',delimiter=' ')
    head,points2 = csv_parse.read('../Baltay-fibers_residual.csv',delimiter=' ')
    points2 = rec.drop_fields(points2,('r','theta'))
    points2['Number'] += 10000 # to distinguish them from the "randoms"
    points = np.hstack((points1,points2))
    return points
#...

def match_positions(data,points):
    """
    Match to the positions in Baltay's files.
    """
    x = 'Y_IMAGE_DBL'
    y = 'X_IMAGE_DBL'

    transform = fvmapping.Transform(points,np.array(zip(data[x],data[y])))
    transform()
    return transform.check()['dist']
#...

def naturalsort(l):
    """Sort the given iterable in the way that humans expect.
    Stolen from:
        http://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python"""
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)
#...

def plot_hist(dist,flux):
    """
    Plot histograms of the distances, labeled appropriately.
    """
    colors = {400:'cyan',4000:'red',10000:'blue',40000:'green',-1:'black'}
    for name in naturalsort(dist.keys()):
        counts,edges = np.histogram(dist[name],range=(0,0.05))
        centers=(edges[1:]+edges[:-1])/2
        if '4_' in name:
            marker = '--'
        else:
            marker = '-'
        try:
            height = int(name[2:])
            label = name.replace('_','\mu m,\;\mathrm{height:}\,')
            label = '$'+label+',\;\mathrm{flux:}'+str(np.round(flux[name]))+'$'
        except ValueError:
            height = -1
            label = 'full noise'
        pylab.plot(centers,counts,marker,lw=2,color=colors[height],label=label)
    pylab.legend(loc='upper right')
    pylab.plot((1./30,1./30),(0,50),'--',lw=1.5,color='grey')
    pylab.ylim((0,40))
    pylab.xlabel('separation (pixels)')
    pylab.ylabel('Number')
    #pylab.savefig('testimage_matches.png')
#...

def main(argv=None):
    if argv is None: argv = sys.argv[1:]
    from optparse import OptionParser, OptionGroup

    usage = '%prog [OPTIONS] FITSFILES'
    usage += '\nRun Sextractor on files and compare with known positions.'
    parser = OptionParser(usage)
    parser.add_option('-p','--plot',dest='plot',action='store_true',default=False,
                      help='Make some histograms of the results (%default).')
    (opts,args) = parser.parse_args(argv)

    if len(args) == 0:
        parser.error('Need to pass some fits image files!')

    points = read_positions()

    dist = {}
    flux = {}
    for f in args:
        if 'testname' in f:
            name = os.path.splitext(f)[0].split('testimage_')[-1]
        else:
            name = os.path.splitext(os.path.basename(f))[0].split('S')[-1]
        print 'Processing:',name
        measured,filename = run_sex(f)
        make_regions(measured,filename)
        #dist[name] = match_positions(measured,points)
        #print 'Matches found:',len(dist[name])
        flux[name] = np.mean(measured['FLUX_ISO'])        
        print 'N, mean(flux) and sigma(flux):',len(measured),flux[name],np.std(measured['FLUX_ISO'])
        #print

    if opts.plot:
        plot_hist(dist,flux)
#...

if __name__ == "__main__":
    import sys
    sys.exit(main())
