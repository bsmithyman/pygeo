# Dip analysis and filtering toolkit

# ----------------------------------------------------------------------
# Authors:

# |-----| Brendan Smithyman (2010)
# | BRS | bsmithyman@eos.ubc.ca
# |-----| University of British Columbia

# ----------------------------------------------------------------------
# References:
#
#  Gradient Square Tensor based dip estimation
#  after L. J. van Vliet and P. W. Verbeek (1995)
#        Estimators for orientation and anisotropy in digitized images
#
#  Folded symmetric dip filter B (variant of wavekilll filter) and its inverse
#  after D. Hale (2007)
#        Local dip filtering with directional Laplacians

# ----------------------------------------------------------------------
# Description:
#  This toolkit implements a number of dip-filtering related algorithms from
#  recent image-processing and seismic literature. This was done with the aim
#  of providing stable local-dip preconditioners for full-waveform inversion of #  seismic data.

# ----------------------------------------------------------------------
# Requirements:
#  Uses imported functions from scipy (might be possible to implement with only
#  numpy) toolkit.

# Optional:
#  The test suite uses plotting functions from matplotlib (pylab)

# ----------------------------------------------------------------------
# Version History:

#  BRS: August - September, 2010
#    Initial implementation of van Vliet and Verbeek dip estimation and
#    Hale B, B^-1 filters.

# ======================================================================
# Library Implementation
# ======================================================================

# Module exports
__all__ = ['orientation', 'orient', 'wkf', 'invwkf', 'rorient', 'diptransfer', 'rdiptransfer']

# Global module imports

import numpy as _np
cimport numpy as _np

import scipy as _sp
import scipy.ndimage as _ndimage

_gf = _ndimage.gaussian_filter
_gf1d = _ndimage.gaussian_filter1d

ctypedef F64_t F64_t

# ----------------------------------------------------------------------
cdef _gsi (_np.ndarray[F64_t, ndim=2] im):
  '''
  Finds the orientation of features (implements GST dip estimation from
  van Vliet and Verbeek (1995). This determines the local dip based on the
  eigenvalue decomposition of the 2D gradient squared tensor.
  '''
  cdef _np.ndarray[F64_t, ndim=2] grads0, grads1

  grads0 = _gf1d(im, 1, 1, 1)
  grads1 = _gf1d(im, 1, 0, 1)
  return [grads1**2, grads0*grads1, grads0**2]

# ----------------------------------------------------------------------
cpdef orientation (_np.ndarray[F64_t, ndim=2] im, gfsize=1):
  '''
  orientation(im, gfsize=1) -> [[lambda1, lambda2],[ori1,ori2]]

  Finds the orientation of features (implements GST dip estimation from
  van Vliet and Verbeek (1995). This determines the local dip based on the
  eigenvalue decomposition of the 2D gradient squared tensor.
  '''

  cdef _np.ndarray[F64_t, ndim=2] g0sq, gmul, g1sq, item, gradsamp, gradsdiff

  # Gets GST components from above
  [g0sq,gmul,g1sq] = [_gf(item, gfsize) for item in _gsi(im)]

  # Computes eigenvalues of the tensor (quadratic equation)
  gradsamp = (g0sq + g1sq)/2.
  gradsdiff = _np.sqrt((g0sq-g1sq)**2 + 4*(gmul**2))/2.

  # Eigenvalues
  lambdas = [_np.real(gradsamp+gradsdiff),_np.real(gradsamp-gradsdiff)]

  # Orientations
  dirs = [_np.real(_np.arctan((lambdas[0]-g0sq)/gmul)),_np.real(_np.arctan(gmul/(lambdas[1]-g1sq)))]

  return [lambdas, dirs]

# ----------------------------------------------------------------------
cpdef orient (im, gfsize=1):
  '''
  orient(im, gfsize=1) -> [ori,aniso]

  Wraps orientation(...) to provide a single orientation and the anisotropy
  representation from vV+V. Nan values may be returned in unstable regions.
  '''

  cdef _np.ndarray[F64_t, ndim=2] l1, l2, a1, a2, aniso

  [[l1,l2],[a1,a2]] = orientation(im, gfsize)
  aniso = 1 - (l2/l1)
  return [a1,aniso]

# ----------------------------------------------------------------------
cdef _mnp (imor):
  '''
  Helper function (defines m and p terms for Hale (2007) B and B^-1
  '''

  # Normal components relative to local orientation
  u2 = _np.sin(imor)
  u1 = _np.cos(imor)

  m = (u1 - u2)/2
  p = (u1 + u2)/2

  return [m,p]

# ----------------------------------------------------------------------
cpdef wkf (im, imor):
  '''
  wkf(im, imor) -> 2D Array

  Implementation of B (functionally) from Hale (2007)
  This is a 'folded' version of a symmetric dip filter (matrix square of
  Claerbout's wavekill filter). Designed by Dave Hale to be similar to
  A^TA, a 3x3 symmetric descendent of a 2x2 wavekill filter. This filter (B)
  folds A^TA about its vertical axis to improve its invertability.
  '''

  # Array of m and p values (image dimensions)
  [mar,par] = _mnp(imor)

  # Output image is initialized at the same size as the original
  imout = _np.zeros(im.shape)

  # Apply the finite-difference stencil explicitly (might not be very fast)
  for j in xrange(1,im.shape[0]-1):
    for i in xrange(1,im.shape[1]):
      m = mar[j,i]
      p = par[j,i]

      # Formed per Hale (2007)
      stencil = 2*m*_np.array([[-m,p],[-p,m],[0,0]]) + 2*p*_np.array([[0,0],[-m,p],[-p,m]])
      imout[j,i] = _np.sum(stencil * im[j-1:j+2,i-1:i+1])

  return imout

# ----------------------------------------------------------------------
cdef _tridiag (
	_np.ndarray[F64_t, ndim=1] a,
	_np.ndarray[F64_t, ndim=1] b,
	_np.ndarray[F64_t, ndim=1] c,
	_np.ndarray[F64_t, ndim=1] d,
	_np.ndarray[F64_t, ndim=1] x):
  '''
  Implementation of TDMA (tridiagonal matrix algorithm) aka Thomas algorithm
  '''

  cdef double *A = <double *> a.data
  cdef double *B = <double *> b.data
  cdef double *C = <double *> c.data
  cdef double *D = <double *> d.data
  cdef double *X = <double *> x.data

  cdef Py_ssize_t p, j
  cdef double temp

  # Setting up length and results array
  n = len(b)

  # Initialization (possible divide by zero)
  C[0] = C[0] / B[0]
  D[0] = D[0] / B[0]

  # Elimination
  for j in xrange(1, n):
    temp = 1 / (B[j] - C[j-1] * A[j])
    C[j] = C[j] * temp
    D[j] = (D[j] - D[j-1] * A[j]) * temp

  # Back-substitution
  X[n-1] = D[n-1]
  for j in xrange(n-2,-1,-1):
    X[j] = D[j] - C[j] * X[j+1]

# ----------------------------------------------------------------------
cpdef invwkf (im, imor):
  '''
  invwkf(im, imor) -> 2D Array

  Implementation of B^-1 (functionally) from Hale (2007)
  This is the inverse of the filter (B) wkf(...). This progresses from left to
  right in the image space with the boundary condition that the results array
  outside the image is zero. This uses the tridiagonal matrix solver method
  suggested by Hale (2007), which is simpler than iteratively solving for the
  inverse of A^TA, and requires no iterations.
  '''

  cdef _np.ndarray[F64_t, ndim=2] m, p, imout
  cdef _np.ndarray[F64_t, ndim=1] a, b, c, d, x
  cdef Py_ssize_t n, n2, i, j

  [m,p] = _mnp(imor)

  # The output image is initialized to zeros
  imout = im.copy()

  # For convenience and brevity
  n = im.shape[0]
  n2 = im.shape[1]
  
  # Set up tridiagonal matrix system and solve for i = 0
  # This is done separately to account for the imout[:,-1] = 0
  # boundary condition and avoid conditionals in the later for loop

  # Setting up tridiagonal matrix terms a,b,c and data vector d
  a = _np.zeros(n)
  b = 1.*_np.ones(n)

  # Note that a and c are identical, but may be changed in the TDMA code.
  # Therefore it is *critical* that c is a copy of a, not a reference to a
  for j in xrange(0,n):
    a[j] = 2.*m[j,0]*p[j,0]

  c = a.copy()

  # Boundary conditions
  a[0] = 0.
  c[n-1] = 0.

  # D is initialized as a copy of the image for the first (i = 0) column.
  # This is because of the zero boundary condition on i = -1.
  d = im[:,0].copy()

  # Set up result vector x
  x = _np.zeros(n)

  # Solve for the first column of the image
  _tridiag(a,b,c,d,x)
  imout[:,0] = x

  # Set up tridiagonal matrix system and solve for 0 < i < n columns
  for i in xrange(1, n2):

    # Setting up tridiagonal matrix terms a,b,c and data vector d
    a = _np.zeros(n)
    b = 1.*_np.ones(n)

    # See above note re: a vs. c
    for j in xrange(0,n):
      a[j] = 2.*m[j,i]*p[j,i]
 
    c = a.copy()

    # Data vector is initialized to zeros
    d = _np.zeros(n)

    # Boundary conditions
    a[0] = 0.
    c[n-1] = 0.

    # Initialize data vector for the top and bottom of the column.
    # Implements boundary condition: imout[-1,:] = imout[n,:] = 0
    d[0] = im[0,i] + 4*m[0,i]*p[0,i]*imout[0,i-1] + 2*p[0,i]**2*imout[1,i-1]
    d[n-1] = im[n-1,i] + 2*m[n-1,i]**2*imout[n-2,i-1] + 4*m[n-1,i]*p[n-1,i]*imout[n-1,i-1]
    
    # Initialize data vector for the remainder of the column (0 < j < n)
    for j in xrange(1, n-1):
      d[j] = im[j,i] + 2*m[j,i]**2*imout[j-1,i-1] + 4*m[j,i]*p[j,i]*imout[j,i-1] + 2*p[j,i]**2*imout[j+1,i-1]

    # Solve for each column of the image
    _tridiag(a,b,c,d,x)
    imout[:,i] = x

  return imout

# ----------------------------------------------------------------------
cpdef rorient (im, gfsize=1):
  '''
  rorient(im, gfsize=1) -> [rori,aniso]

  Wraps orientation(...) to provide a single orientation and the anisotropy
  representation from vV+V. Undefined values are converted to random dips
  to avoid returning NaNs.
  '''
  
  [ori, aniso] = orient(im, gfsize)
  aniso = _np.nan_to_num(aniso)
  rpad = _sp.rand(*im.shape) * _np.pi - _np.pi/2
  rori = _np.nan_to_num((1 - _np.isnan(ori))*ori) + _np.isnan(ori) * rpad

  return [rori, aniso]

# ----------------------------------------------------------------------
cpdef diptransfer (source, dest, gfsize=1):
  '''
  diptransfer(source, dest, gfsize=1) -> 2D Array

  Automated dip-transfer routine using orient(...)
  May be unstable in smooth/isotropic images.
  '''

  [ori,aniso] = orient(source, gfsize)
  newimg = invwkf(dest, ori)
  return newimg

# ----------------------------------------------------------------------
cpdef rdiptransfer (source, dest, gfsize=1):
  '''
  rdiptransfer(source, dest, gfsize=1) -> 2D Array

  Automated dip-transfer routine using rorient(...) and aware of
  anisotropy. Attempts to reconstruct a stable dip field using random
  in areas where dip-scanning fails and returning the original image
  in isotropic regions.
  '''

  [ori,aniso] = rorient(source, gfsize)
  [sori,saniso] = rorient(dest, gfsize)
  nori = aniso*ori + (1-aniso)*sori
  newimg = aniso*invwkf(dest, nori) + (1-aniso)*dest
  return newimg
