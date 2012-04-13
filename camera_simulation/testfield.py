"""
Generate four example plots with 2d gaussians at the currently accepted
positions of the 'residual' fibers.

For testing the simCCD code.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import simccd

w = 2*np.sqrt(2*np.log(2)) # FWHM = w*sigma
scale = 6. # microns/pixel


def makefile(simCCD,width,height,skew=False,skew_a=1./5000.):
    '''width in microns.'''
    simCCD(100,layout='Baltay_default',width=width/scale,height=height,skew=skew,skew_a=skew_a)
    filename = 'testimage_'+str(width)+'_'+str(int(np.round(height)))+'.fits'
    simCCD.save(filename,clobber=True)
    print 'saved:',filename
#...

simCCD = simccd.SimCCD(Ndim=(1024,1024))
makefile(simCCD,4,4e2)
makefile(simCCD,4,4e3)
makefile(simCCD,4,1e4)
makefile(simCCD,4,4e4)

makefile(simCCD,5,4e2)
makefile(simCCD,5,4e3)
makefile(simCCD,5,1e4)
makefile(simCCD,5,4e4)

# comment out to make the below plots appear.
sys.exit()


plt.matshow(simCCD.image,cmap=cm.gist_ncar)
plt.title('float image, linear scaling')

plt.matshow(np.log10(simCCD.image),cmap=cm.gist_ncar)
plt.title('float image, log10 scaling')

plt.matshow(simCCD.image_int,cmap=cm.gist_ncar)
plt.title('int image, linear scaling')

plt.matshow(np.log10(simCCD.image_int),cmap=cm.gist_ncar)
plt.title('int image, log10 scaling')
plt.show()

sys.exit()
