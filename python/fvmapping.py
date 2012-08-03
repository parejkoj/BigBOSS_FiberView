"""
Transform points in the image plane to points in the fiber plane.
"""
import numpy as np

from kdtree_test.naive_pairs import xymatch

class Transform:
    """
    Transform a point in the xy image plane, xi, to its equivalent point
    in the fiber plane, xf.
    
    This proces includes finding matching points, creating a first
    guess at the transformation, and then minimizing the full point fit.
    """
    def __init__(self,iPoints,fPoints):
        """
        Provide points to find the transform T from i to f: T*i -> f
        """
        if iPoints.shape != fPoints.shape:
            raise ValueError('Both sets of coordinates must have the same size!')
        self.iPoints = iPoints
        self.fPoints = fPoints
        self.T = np.identity(3) # identity for an affine transformation
    #...

    def __call__(self):
        """Find the best fitting transformation that takes iPoints to fPoints."""
        self._initial_fit()
        match = self.check()
        distmean = match['dist'].mean()
        if distmean > 10:
            raise ValueError('Initial guess failed: mean transform separation is %4f'%distmean)
        self._fit(self.iPoints[match['idx1']],self.fPoints[match['idx2']])
    #...

    def check(self):
        """
        Check how far apart the transformed points are from the real points.
        Returns a matching array: (idx_i,idx_f,dist)
        """
        points = self.transform_all()
        match = xymatch(points,self.fPoints,100.,maxmatch=1)
        return match
    #...
    
    def transform_all(self):
        """Return the transformed locations of all image plane points."""
        points = np.empty(self.iPoints.shape)
        for i,pi in enumerate(self.iPoints):
            points[i] = self.transform(pi)
        return points
    #...

    def transform(self,xi):
        """
        Transform a point in the xy image plane, xi, to it's equivalent point
        in the fiber plane, xf.
        """
        # NOTE: remember that T is an affine 3x3 matrix. So we need to convert xi
        # to a 3x1 "affine vector" and chop off the "1" in the last row before returning.
        return np.dot(self.T,self._affine_vector(xi))[:2]
    #..

    def get_scale(self):
        """Return the overall scale factor of the matrix."""
        x0 = np.array([1.,1.])
        x0 /= np.linalg.norm(x0)
        x1 = self.transform(x0)
        return np.linalg.norm(x1)
    #...

    def _affine_vector(self,x):
        """Returns the 3x1 affine transform vector (x1,x2,1)."""
        return np.array((x[0],x[1],1.))
    #...
    
    def _make_translation(self,vector):
        """Return a translation matrix for the given vector."""
        C = np.identity(3)
        C[:2,2] = vector
        return C
    #...
    
    def _make_linear_transform(self,M):
        """Return a 3x3 linear transformation matrix, given a 2x2 matrix."""
        Mnew = np.identity(3)
        Mnew[:2,:2] = M
        return Mnew
    #...

    def _initial_fit(self):
        """
        Initial guess at the transformation: xf = T*xi using the
        furthest pair of points.
        """
        self._furthest()
        self._fit(self.iPoints[self.iFar],self.fPoints[self.fFar])
    #...

    def _fit(self,iP,fP):
        """
        Compute the best-fitting transformation: xf = T*xi using SVD.
        iP and fP are the points to use when computing the transformation,
        and are assumed to be in matched order:
            T(iP[0]) -> fP[0],
            T(iP[1]) -> fP[0], etc...

        Here is a mathematical summary of the steps:
            igl.ethz.ch/projects/ARAP/svd_rot.pdf
        """
        # Find the centroids of the data
        iCent = iP.mean(axis=0)
        fCent = fP.mean(axis=0)
        # shift the data to the origin
        iNew = iP - iCent
        fNew = fP - fCent

        # NOTE: the below flips the image around the y-axis to correct for
        # the known flip in the image coordinates. This doesn't appear to be
        # picked up by the SVD.
        #yflip = np.identity(2)
        #yflip[1,1] = -1
        #iNew = np.dot(yflip,iNew.T).T

        # find the optimal rotation+scaling via SVD
        H = np.zeros((2,2))
        for i in range(iP.shape[0]):
            H += np.outer(iNew[i,:],fNew[i,:])
        # remember that numpy's SVD gives H=U*S*V, not U*S*V.T, as the,
        # e.g., wikipedia definition. So we need to take V.T when using it in
        # the matrix product below.
        U,S,V = np.linalg.svd(H)
        M = np.dot(V.T,U.T)
        det = np.identity(2)
        det[1,1] = np.linalg.det(M)
        M = np.dot(np.dot(V.T,det),U.T)
        # create 3x3 affine transform matrix.
        Ci = self._make_translation(-iCent)
        Cf = self._make_translation(fCent)
        Mnew = self._make_linear_transform(M)
        self.T = np.dot(np.dot(Cf,Mnew),Ci)
        #import pdb
        #pdb.set_trace()
    #...

    def _minimize(self):
        """Minimize the difference between the transformed fiber positions and the actual image."""
        return 0
    #...

    def _furthest(self):
        """Find the two points separated by the largest distance in the input lists."""

        def find_furthest(points):
            """Find the furthest pair in the array points."""
            r = 0 # r^2
            for i in range(len(points)):
                for j in range(i,len(points)):
                    ri = ((points[i]-points[j])**2).sum()
                    if ri > r:
                        far = [i,j]
                        r = ri
            return np.sqrt(r),far
        #...
        self.iR,self.iFar = find_furthest(self.iPoints)
        self.fR,self.fFar = find_furthest(self.fPoints)
    #...
#...
