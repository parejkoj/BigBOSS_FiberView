"""
Generate four example plots with 2d gaussians at the currently accepted
positions of the 'residual' fibers.

For testing the simCCD code.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
# import simCCD
import simccd2

# inside a python window, do this:

simCCD = simccd.SimCCD(Ndim=(1000,1000))
simCCD(100,layout='Baltay_default')

imcopy = simCCD.image.copy() # save a copy of the image
simCCD.dark_current()

# see how it looks. (figure out how to best to do this!)
# if it didn't work, replace the image and try again.

simCCD.image = imcopy
simCCD.dark_current()

# for checking the execution time of a function, look at the timeit function
# also, look at stackexchange for help with optimizing numpy/python
# and look at the links I emailed you for ideas.
# I bet you can get a ~5x speed up in shot_noise without much effort.

# do some statistical tests to see if your added noise is appropriate.

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

# for i in range(1000):
    # for j in range(1000):
        # image[i][j] +=np.random.normal(12.5)
        #what is the standard deviation for readnoise?
# image *= 12.5

for c in coords:
    image += gaussian_profile(2, c[0], c[1],3,xygrid=xygrid)
    
# def add_darkcurrent(x,temperature, duration):
    # for i in range(1000):
        # for j in range(1000):
            # if temperature is 296:
                # x[i][j] += 15*duration
            # else:
                # print "What is the temperature?"
                
def add_darkcurrent_roomtemp(x, duration):
    for i in range(1000):
        for j in range(1000):
            x[i][j] += 15*duration
                
add_darkcurrent(image,10)

for i in range(1000):
    for j in range(1000):
        image[i][j] += np.random.normal(12.5)


    # image += np.random.poisson(3,gaussian_profile(2, c[0], c[1],xygrid=xygrid))
# print image
k
