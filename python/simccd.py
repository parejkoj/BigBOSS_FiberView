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
import pyfits
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
    def __init__(self,sigma,skew=False,height=40000,skew_a=None): #skew_param=0,skew_scale=0,
        self.sigma = sigma
        # TBD: what is the correct prefactor here?
        self.height = height#*(sigma*sqrt(2*pi))**2
        self.skew = skew
        if skew_a is None:
            self.skew_a = 1/(sqrt(2)*500)
        else:
            self.skew_a = skew_a
        # self.skew_param = skew_param
        # self.skew_scale = skew_scale
    #...

    def skewness(self,center_x,center_y):
        """Returns a radially-dependent skewness parameter, alpha."""
        return self.skew_a*sqrt(center_x**2+center_y**2)

    def cum_dist_func(self,x,y,center_x,center_y,skew):
        return (0.5*sqrt(pi)*self.sigma)**2*(1-erf(skew*(center_x-x)/self.sigma))*((1-erf(skew*(center_y-y)/self.sigma)))

    def gauss(self,x,y,center_x,center_y):
        """Returns a symmetric 2d gaussian profile, evaluated on the grid x,y."""
        # TBD: check that the prefactor is correct.
        return self.height*1./(self.sigma*sqrt(2*pi))*np.exp(-(((center_x-x)/self.sigma)**2+((center_y-y)/self.sigma)**2)/2)

    def skew_gaussian(self,x,y,center_x,center_y):
        """Returns a radially skewed 2d gaussian, evaluated on the grid x,y."""
        skew = self.skewness(center_x,center_y)
        return (2./self.sigma)*self.gauss(x,y,center_x,center_y)*self.cum_dist_func(x,y,center_x,center_y,skew)

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

        if self.skew is True:
            return self.skew_gaussian(x,y,center_x,center_y)
        else:
            return self.gauss(x,y,center_x,center_y)
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

    def __call__(self,Npoints=None,layout='Baltay_default',resid_only=False,
                 width=2.,height=40000.,skew=False,skew_a=None):
        """
        Generates an image with Npoints distributed by layout.
        Valid options for layout are:
            xygrid
            random
            Baltay_default
            ???
        set resid_only to True to only include the 'residual' points in the output.
        width is the width of the Gaussian (==RMS), in pixels.
        height is the value of the peak of the Gaussian.
        """
        self.image[:] = 0
        self.image_int[:] = 0
        gaussian = Gaussian(width,height=height,skew=skew,skew_a=skew_a)
        if layout == 'Baltay_default':
            coordinates = np.loadtxt('../Baltay-fibers_residual.csv')
            for c in coordinates:
                result = gaussian(c[3],c[4],xygrid=self.xygrid)
                self.image += result
            if not resid_only:
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
    
    def save(self,filename,clobber=False):
        """
        Save the current image to a .fits file named filename.

        Set clobber=True to overwrite the file if it exists."""
        hdu = pyfits.PrimaryHDU(self.image_int)
        hdu.writeto(filename,clobber=clobber)
    #...
#...
