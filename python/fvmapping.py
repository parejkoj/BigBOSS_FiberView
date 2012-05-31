"""
Transform points in the image plane to points in the fiber plane.
"""
import numpy as np

class Transform:
    """
    Transform a point in the xy image plane, xi, to it's equivalent point
    in the fiber plane, xf.
    
    This proces includes finding matching residual points, creating a first
    guess at the transformation, and then minimizing the residual point fit.
    """
    def __init__(self,iResid,fResid):
        """
        Takes the coordinates of the 'residuals' in each plane.
        """
        self.iResid = iResid
        self.fResid = fResid
        self.A = np.identity(2) # identity for scale/rotation
        self.b = 0 # no translation
    #...

    def __call__(self):
        """Find the best fitting transformation that takes iResid to fResid."""
        self._furthest()
        self._initial_fit()
        self._minimize()
    #...

    def _initial_fit(self):
        """Generate an initial guess at the fit using the furthest pair of points."""
        _transform
    #...

    def _minimize(self):
        """Minimize the difference between the transformed fiber positions and the actual image."""
        return 0
    #...

    def _furthest(self):
        """Find the two points separated by the largest distance in the input lists."""
        # All points in the image plane are > 0
        # so just search for the two with largest and smallest radius
        r2 = (self.iResid**2).sum(axis=1)
        self.iFar = (r2 == r2.min()) | (r2 == r2.max())
        self.iR = np.sqrt((self.iResid[self.iFar]**2).sum())

        # the fiber plane can have points in all 4 quadrants
        # so we have to search through all the pairs
        self.fR = 0
        for i in range(len(self.fResid)):
            for j in range(i,len(self.fResid)):
                ri = ((self.fResid[i]-self.fResid[j])**2).sum()
                if ri > self.fR:
                    self.fFar = [i,j]
                    self.fR = ri
        self.fR = np.sqrt(self.fR)
    #...

    def _transform(self,xi):
        """
        Transform a point in the xy image plane, xi, to it's equivalent point
        in the fiber plane, xf.
        """
        return np.dot(self.A,xi)+self.b
    #...
#...
