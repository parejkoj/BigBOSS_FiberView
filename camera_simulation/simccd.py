"""
this is the framework for generating simulated CCD output.

Generates images with points distributed in Gaussian profiles, with errors
like those expected from the FLI PL50100 dark current and read noise.

Example:
    simCCD = SimCCD(Ndim=(1000,1000))
    simCCD(100,layout='Baltay_default')
    pylab.matshow(simCCD.image) # plot the array
    simCCD.save('xy_simulation.fits')
    simCCD(50,layout='random')
    simCCD.save('random_simulation.fits')
"""

# import numpy as np
# from numpy import *
# from pylab import *

# import pyfits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from scipy import linspace
from scipy import pi,sqrt,exp
from scipy.special import erf

# NOTE: scipy.ndimage and scipy.signal may have some useful tools.
# NOTE: for testing purposes, pylab.matshow could be useful.

class Gaussian:
    """Class to make single Gaussian Profile"""
    def __init__(self,sigma,skew_param,skew_scale,height=40000):
        self.sigma = sigma
        self.height = height
        self.skew_param = skew_param
        self.skew_scale = skew_scale
        
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
            return height*np.exp(-(((center_x-x)/sigma)**2+((center_y-y)/sigma)**2)/2)

        def cum_dist_func(x,y,height,sigma):
            return (-0.5*sqrt(pi)*height*sigma)*(1-erf((center_x-x)/sigma))*((1-erf((center_y-y)/sigma))
        
        def skew_gaussian(x,y,w=1,a):
            tx=x/w
            ty=y/w
            def skewness(a)
                return lambda a: 3*sqrt(center_x^2+center_y^2)
                #enter function, parameters which govern radial distortions
            return (2/w)*gauss(tx,ty,self.height,self.sigma)*cum_dist_func(a*tx,a*ty,self.height.self.sigma)
       
        return skew_gaussian(x,y)
            
        # return gauss(x,y,self.height,self.sigma)
        
        
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
        self.image = np.random.poisson(self.image)
    #...

    def read_noise(self):
        """
        Adds read noise to every pixel value, using the parameters
        of the Kodak CCD.
        """
        self.image += np.random.normal(12.5,size=self.image.size).reshape(self.image.shape)
    #...
    
    def dark_current(self,t=1,T=296):
        """
        Adds dark current to every pixel value, assuming temperature T
        (in degrees Celcius) and exposure duration t (in seconds).
        """
        self.image += np.random.poisson(15*t,self.image.size).reshape(self.image.shape)
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
            coordinates = np.loadtxt('../Baltay-fibers_residual.csv')
            for c in coordinates:
                self.image += gaussian(c[3],c[4],xygrid=self.xygrid)
            coordinates = np.loadtxt('../Baltay-fibers_random.csv')
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