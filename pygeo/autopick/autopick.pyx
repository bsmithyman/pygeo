import cython
cimport cython
import numpy as np
cimport numpy as np

from scipy.fftpack import hilbert as _hilbert
from scipy.ndimage import median_filter as _mf

#np.import_array()

ctypedef np.float32_t F32_t

# ------------------------------------------------------------------------
# Configuration Options

eardamp = 100
earwindow = 40

# ------------------------------------------------------------------------
# C External References
cdef extern void c_energyRatio "energyRatio" (F32_t inarr[], F32_t outarr[], Py_ssize_t arrL, Py_ssize_t arrW, Py_ssize_t strideL, Py_ssize_t strideW, Py_ssize_t windowsize, double damp) nogil

def energyRatio (np.ndarray[F32_t, ndim=2] traces, Py_ssize_t windowsize=earwindow, double damp=eardamp):
  '''
  Generates a STA/LTA (Short-Term Average over Long-Term Average) filtered
  result that highlights sharp increases in trace amplitude.
  '''

  cdef np.ndarray[F32_t, ndim=2] result
  result = np.empty((traces.shape[0],traces.shape[1]), dtype=np.float32)

  cdef Py_ssize_t arrL, arrW, strideL, strideW
  cdef F32_t *inarr = <F32_t *> traces.data
  cdef F32_t *outarr = <F32_t *> result.data

  arrL = traces.shape[0]
  arrW = traces.shape[1]
  strideL = traces.strides[0]
  strideW = traces.strides[1]

  #with cython.nogil:
  c_energyRatio(inarr, outarr, arrL, arrW, strideL, strideW, windowsize, damp)

  return result

def envelope (np.ndarray[F32_t, ndim=2] traces):
  '''
  Returns the envelope of a time series.
  '''

  cdef Py_ssize_t i
  cdef np.ndarray[F32_t, ndim=2] result
  result = np.empty((traces.shape[0],traces.shape[1]), dtype=np.float32)

  for i in xrange(traces.shape[0]):
    result[i,:] = np.sqrt(traces[i,:]**2 + _hilbert(traces[i,:])**2)

  return result

def _argsig (traces, count):
  locs = traces.argsort(axis=1)[:,-count:]
  res = np.median(locs, axis=1)
  return res

def xcorval (traces, offset=10, medianlen=10):

  # FFT Functions to use
  fft = np.fft.rfft
  ifft = np.fft.irfft
  fftshift = np.fft.fftshift

  # Create spectra for the original and time-reversed traces
  nT = fft(traces)
  rT = fft(np.fliplr(traces))

  # Compute the crosscorrelation and shift the spike
  offtr = np.zeros(nT.shape, dtype=nT.dtype) + traces.shape[1]/2
  offtr[offset/2:-offset/2] = nT[offset:]*rT[:-offset]
  shifttr = fftshift(ifft(offtr), axes=1)

  # Compute the envelope of the xcf and choose the centre
  shifts = _argsig(envelope(shifttr.astype(np.float32)), medianlen) - traces.shape[1]/2.
  
  return _mf(shifts / offset, 3), shifttr

def xcorvalize (traces, offrange, medianlen=10):
  x = np.zeros(len(traces))
  noff = offrange[1] - offrange[0]

  for offset in xrange(*offrange):
    x += xcorval(traces, offset, medianlen)

  x = x / noff

  y = np.empty(x.shape)
  z = 0.

  for i in xrange(len(x)):
    z += x[i]
    y[i] = z

  return y
