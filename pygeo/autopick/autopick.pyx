import numpy as np
cimport numpy as np
import cython
cimport cython

np.import_array()

ctypedef np.float32_t F32_t

# ------------------------------------------------------------------------
# Configuration Options

eardamp = 100
earwindow = 40

# ------------------------------------------------------------------------
# C External References
cdef extern void c_energyRatio "energyRatio" (F32_t inarr[], F32_t outarr[], Py_ssize_t arrL, Py_ssize_t arrW, Py_ssize_t strideL, Py_ssize_t strideW, Py_ssize_t windowsize, double damp) nogil

def energyRatio (np.ndarray[F32_t, ndim=2] traces, Py_ssize_t windowsize=earwindow, double damp=eardamp):

  cdef np.ndarray[F32_t, ndim=2] result
  result = np.empty((traces.shape[0],traces.shape[1]), dtype=np.float32)

  cdef Py_ssize_t arrL, arrW, strideL, strideW
  cdef F32_t *inarr = <F32_t *> traces.data
  cdef F32_t *outarr = <F32_t *> result.data

  arrL = traces.shape[0]
  arrW = traces.shape[1]
  strideL = traces.strides[0]
  strideW = traces.strides[1]

  with cython.nogil:
    c_energyRatio(inarr, outarr, arrL, arrW, strideL, strideW, windowsize, damp)

  return result

