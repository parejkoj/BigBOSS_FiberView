#/usr/bin/env python
"""
this is the framework for generating simulated CCD output.

Generates images with points distributed in Gaussian profiles, with errors
like those expected from the FLI PL50100 dark current and read noise.

Example:
    simCCD = SimCCD(Ndim=(1000,1000))
    simCCD(100,layout='xygrid')
    pylab.matshow(simCCD.image) # plot the array
    simCCD.save('xy_simulation.fits')
    simCCD(50,layout='random')
    simCCD.save('random_simulation.fits')
"""

# import numpy as np
# from numpy import *
# from pylab import *

import pyfits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# NOTE: scipy.ndimage and scipy.signal may have some useful tools.
# NOTE: for testing purposes, pylab.matshow could be useful.

class Gaussian:
    """Class to make single Gaussian Profile"""
    def __init__(self,sigma,height=40000):
        self.sigma = sigma
        self.height = height
        
    def __call__(self,center_x,center_y,xygrid=None):
        """
        Returns a symmetric 2d Gaussian profile, with a sigma in pixels.
        Centered in the middle of a gridsize x gridsize array.
        The height is the value at the peak of the Gaussian.
        """
        if xygrid is None:
            x,y = np.mgrid[0:20,0:20]
        else:
            x,y = xygrid

        def gauss(x,y,height,sigma):
            return height*np.exp(-(((center_x-x)/self.sigma)**2+((center_y-y)/self.sigma)**2)/2)

        return gauss(x,y,self.height,self.sigma)
    #...
#...

class SimCCD:
    """Generate a simulated CCD image."""
    def __init__(self,Ndim=(1024,1024),dark=0,readnoise=0):#,shot_time=1): 
        """
        Ndim a tuple of the dimensions of the resulting image.
        """
        self.Nx = Ndim[0]
        self.Ny = Ndim[1]
        self.image = np.zeros([self.Nx,self.Ny],dtype=np.float64)
        self.image_int = np.empty(self.image.shape,dtype=np.uint16)
        self.dark = dark
        self.readnoise = readnoise
        # self.shot_time = shot_time
        x,y = np.mgrid[0:self.Nx,0:self.Ny]
        x -= self.Nx/2
        y -= self.Ny/2
        self.xygrid = x,y
    #...    

    def shot_noise(self):
        """
        Replace every pixel value in self.image with a Poisson random
        variable drawn from its own distribution.
        """
        self.image[i][j] = [map(np.random.poisson,range(self.Nx))][map(np.random.poisson,range(self.Ny))]
        
        # for i in range(self.Nx):
            # for j in range(self.Ny):
                # self.image[i][j] = np.random.poisson(self.image[i][j])
    #...

    def read_noise(self):
        """
        Adds read noise to every pixel value, using the parameters
        of the Kodak CCD.
        """
        self.image[i][j] = [map(+ np.random.normal(12.5),range(self.Nx))][map(+ np.random.normal(12.5),range(self.Ny))]
                # self.image[i][j] = [map(+= np.random.normal(12.5),range(self.Nx))][map(+= np.random.normal(12.5),range(self.Ny))]

        # for i in range(self.Nx):
            # for j in range(self.Ny):
                # self.image[i][j] += np.random.normal(12.5)
    #...
    
    def dark_current(self,t=1,T=296):
        """
        Adds dark current to every pixel value, assuming temperature T
        (in degrees Celcius) and exposure duration t (in seconds).
        """        
        self.image [i][j] = [map(+ 15*t,range(self.Nx))][map(+ 15*t,range(self.Ny))]
        # self.image [i][j] = [map(+= 15*t,range(self.Nx))][map(+= 15*t,range(self.Ny))]

        for i in range(self.Nx):
            for j in range(self.Ny):
                self.image[i][j] += 15*t
    #...

    def __call__(self,Npoints=None,layout='Baltay_default'):
        """
        Generates an image with Npoints distributed by layout.
        Valid options for layout are:
            xygrid
            random
            Baltay_default
            ???
        """
        gaussian = Gaussian(2)
        if layout == 'Baltay_default':
            coordinates = np.loadtxt('Baltay-fibers_residual.csv')
            # coordinates = [[0.000,0.000],
         # [0.000,-120.000],
         # [103.923,-60.000],
         # [103.923,60.000],
         # [0.000,120.000],
         # [207.846,-120.000],
         # [240.000,0.000],
         # [207.846,120.000],
         # [120.000,207.846],
         # [0.000,240.000],
         # [354.531,-62.513],
         # [354.531,62.514],
         # [311.769,180.000],
         # [231.403,275.776],
         # [463.644,-124.233],
         # [480.000,0.000],
         # [463.644,124.233],
         [415.692,240.000]]
            # np.loadtxt('Baltay-fibers_residual.csv')
            for c in coordinates:
                self.image += gaussian(c[3],c[4],xygrid=self.xygrid)
            coordinates = np.loadtxt('Baltay-fibers_random.csv')
        # coordinates = [[427.580,26.662],
        # [38.179,139.204],
        # [113.216,-19.674],
        # [80.793,-0.737],
        # [168.867,-45.043],
        # [24.339,194.961],
        # [291.616,-16.767],
        # [53.059,54.082],
        # [284.854,-85.151],
        # [119.787,-96.488],
        # [358.881,92.067],
        # [479.870,27.114],
        # [154.306,-84.022],
        # [207.077,237.172],
        # [337.906,234.332],
        # [393.574,261.788],
        # [426.812,5.605],
        # [409.772,0.227],
        # [95.702,35.513],
        # [182.104,4.770],
        # [190.652,-110.737],
        # [134.446,107.362],
        # [238.088,-81.042],
        # [123.950,32.216],
        # [165.289,218.759],
        # [67.968,188.407],
        # [179.692,79.141],
        # [154.166,162.511],
        # [212.239,169.032],
        # [371.963,263.334],
        # [154.042,-46.999],
        # [123.555,150.778],
        # [192.932,20.614],
        # [154.550,149.070],
        # [76.317,29.809],
        # [352.885,43.599],
        # [417.108,-44.411],
        # [113.831,253.659],
        # [86.202,-107.610],
        # [336.788,-39.042],
        # [9.031,207.260],
        # [101.543,-37.310],
        # [241.931,38.295],
        # [348.878,101.515],
        # [454.260,103.724],
        # [45.129,197.040]]
            # np.loadtxt('Baltay-fibers_random.csv')
            for c in coordinates:
                self.image += gaussian(c[1],c[2],xygrid=self.xygrid)
        elif layout == 'xygrid':
            raise ValueError("xygrid is not yet implemented. Fix this!")
        else:
            raise ValueError("I don't understand this layout: "+layout)

        # TBD: check that this order is correct.
        # TBD: do we apply shot noise to pixels that only have dark current?
        self.dark_current()
        self.shot_noise()
        self.read_noise()

        # The CCD output will be 16 bit unsigned integers.
        np.round(self.image,out=self.image_int)
    #...
    
    def save(self,filename):
        """Save the current image to a .fits file named filename."""
        hdu = pyfits.PrimaryHDU(self.image_int)
        hdu.writeto(filename)
    #...
#...
