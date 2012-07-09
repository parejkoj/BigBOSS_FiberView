#!/usr/bin/env python
"""
Test the fvmapping code.
"""
import numpy as np
import pyfits
import csv_parse

import fvmapping
"""
# TBD: eventually could make this read several files, or user option
image = pyfits.open('../data/residual_only_4_1e4_sex.fits')[1].data
image = np.array(zip(image['X_IMAGE_DBL'],image['Y_IMAGE_DBL']))
head,fibers = csv_parse.read('../Baltay-fibers_residual.csv',delimiter=' ')
fibers = np.array(zip(fibers['x'],fibers['y']))

transform = fvmapping.Transform(image,fibers)
transform()
print transform.check()
"""
t = np.radians(45.) # 30 degrees
R = np.array(((np.cos(t),-np.sin(t)),(np.sin(t),np.cos(t))))
S = np.identity(2)
S[0] *= 2
S[1] *= 3
print S
print R
M = np.dot(R,S)
M = R
xi = np.array(((0,0),(1,1),(2,3),(4,9)),dtype='f8')
xf = np.array([np.dot(M,x) for x in xi]) + np.array((1,2))
transform = fvmapping.Transform(xi,xf)
transform()
#print transform.check()
# pure translation
for i in range(0,40):
    print i,
    xi = np.array(((0,0),(1,1),(2,3),(4,i)),dtype='f8')
    xf = np.array([np.dot(M,x) for x in xi]) + np.array((1,2))
    transform = fvmapping.Transform(xi,xf)
    transform()
    print transform.check()
