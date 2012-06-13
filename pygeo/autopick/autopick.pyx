
# pygeo - a distribution of tools for managing geophysical data
# Copyright (C) 2011, 2012 Brendan Smithyman

# This file is part of pygeo.

# pygeo is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# pygeo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with pygeo.  If not, see <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------

import cython
cimport cython
import numpy as np
cimport numpy as np

from scipy.fftpack import hilbert as _hilbert
from scipy.ndimage import median_filter as _mf

#np.import_array()
ctypedef np.float32_t F32_t

import pyopencl as cl
import struct

clfilename = 'filter.cl'


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

def energyRatioCL (np.ndarray[F32_t, ndim=2] traces, Py_ssize_t windowsize=earwindow, double damp=eardamp):
  '''
  Generates a STA/LTA (Short-Term Average over Long-Term Average) filtered
  result that highlights sharp increases in trace amplitude.
  '''

  cdef np.ndarray[F32_t, ndim=2] result
  result = np.empty((traces.shape[0],traces.shape[1]), dtype=np.float32)

  cdef Py_ssize_t arrL, arrW, strideL, strideW, jsize, global_size, local_size

  arrL = traces.shape[0]
  arrW = traces.shape[1]
  strideL = traces.strides[0]
  strideW = traces.strides[1]
  jsize = strideL/4

  ctx = cl.create_some_context()
  queue = cl.CommandQueue(ctx)
  mf = cl.mem_flags

  with open(clfilename, 'r') as fp:
    fstr = ''.join(fp.readlines())
  program = cl.Program(ctx, fstr).build()

  traces_buf = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=traces)
  result_buf = cl.Buffer(ctx, mf.WRITE_ONLY, result.nbytes)

  params = struct.pack('iiif', arrW, jsize, windowsize, damp)
  params_buf = cl.Buffer(ctx, mf.READ_ONLY, len(params))
  cl.enqueue_write_buffer(queue, params_buf, params).wait()

  program.energyRatio(queue, (arrL,), None, traces_buf, result_buf, params_buf)
  queue.finish()
  
  cl.enqueue_read_buffer(queue, result_buf, result).wait()

  return result

