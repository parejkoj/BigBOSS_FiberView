#!/usr/bin/env python
"""
Test the fvmapping code.
"""
import pyfits
import csv_parse

import fvmapping

# TBD: eventually could make this read several files, or user option
image = pyfits.open('../data/residual_only_4_1e4_sex.fits')[1].data
image = np.array(zip(image['X_IMAGE_DBL'],image['Y_IMAGE_DBL']))
head,fibers = csv_parse.read('../Baltay-fibers_residual.csv',delimiter=' ')
fibers = np.array(zip(fibers['x'],fibers['y']))

transform = fvmapping.Transform(image,fibers)
transform()
