#!/usr/bin/env python
"""
Overplot the sextractor regions identified in an image in a ds9 window.
"""
import glob
import os
import pyfits

import ds9

def main(argv=None):
    if argv is None: argv = sys.argv[1:]
    from optparse import OptionParser, OptionGroup

    usage = '%prog [OPTIONS] FILE1.fits [FILE2.fits ...]'
    usage += '\nPlot each file in a ds9 instance, with optional regions overplotted.'
    parser = OptionParser(usage)
    (opts,args) = parser.parse_args(argv)

    if len(args) == 0:
        parser.error('Need to pass some fits image files!')

    d = ds9.ds9()

    # TBD: just one file at a time for now.
    f = args[0]
    image = pyfits.open(f)[0].data
    d.set_np2arr(image)
    d.set('cmap b')
    d.set('scale log')
    d.set('scale limits 2000 20000')

    regfiles = '../data/sex/'+os.path.splitext(os.path.basename(f))[0]+'*.reg'
    regions = glob.glob(regfiles)
    for r in regions:
        d.set('regions load '+r)

    import pdb
    pdb.set_trace()
#...

if __name__ == "__main__":
    import sys
    sys.exit(main())
