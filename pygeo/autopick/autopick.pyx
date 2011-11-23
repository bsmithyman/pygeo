import numpy as _np
cimport numpy as _np
cimport cython

import scipy as _sp
import scipy.ndimage as _ndimage
import scipy.fftpack as _fftpack
import inspect as _inspect
import multiprocessing as _mp

from pygeo.segyread import SEGYFile

ctypedef _np.float32_t F32_t
ctypedef _np.float64_t F64_t

# ------------------------------------------------------------------------
# Configuration Options

eardamp = 100
earwindow = 40

# ------------------------------------------------------------------------
# Functions
def _isDescendant (obj, classname):
  return classname in [item.__name__ for item in _inspect.getmro(obj.__class__)]

def _window (ints, wstart=10, function=_np.mean):
  outs = _np.zeros(ints.shape)
  for i in xrange(wstart,len(ints)):
    outs[i] = function(ints[i-wstart:i])

  return outs

def _windowStart (ints, function=_np.mean):
  outs = _np.zeros(ints.shape)
  for i in xrange(len(ints)):
    outs[i] = function(ints[:i])
  return outs

def _energy (ints):
  return _np.sum(_np.square(ints))

def _energyRatioFilter (options):
  [ints, windowsize, damp] = options
  b = _window(ints,windowsize,_energy)
  B = _windowStart(ints,_energy)
  return b/(B+damp)

@cython.wraparound(False)
@cython.boundscheck(False)
def energyRatio (traces, windowsize=earwindow, damp=eardamp, pool=None):
  '''
  Return the energy ratio between a sliding window of a given size and the
  total energy to that point in the file.
  '''

  if (len(traces.shape) == 1):
    traces = _energyRatioFilter(traces, windowsize, damp)
  else:
    if (pool is None):
      try:
        pool = _mp.Pool()
      except:
        myMap = map
      else:
        myMap = pool.map

    nt = len(traces)

    traces = _np.array(myMap(_energyRatioFilter, zip(traces, nt*[windowsize], nt*[damp])))

  return traces

# ------------------------------------------------------------------------
# Classes
class SEGYPicker (SEGYFile):
  '''
  Descendant of SEGYFile optimized and extended for automatic first-arrival
  picking.
  '''

  # Override Options
  #verbose = False

  def _window (self, ints, wstart=10, function=_np.mean):
    outs = _np.zeros(ints.shape)
    for i in xrange(wstart,len(ints)):
      outs[i] = function(ints[i-wstart:i])

    return outs

  def _windowStart (self, ints, function=_np.mean):
    outs = _np.zeros(ints.shape)
    for i in xrange(len(ints)):
      outs[i] = function(ints[:i])
    return outs

  def _energy (self, ints):
    return _np.sum(_np.square(ints))

  def _energyRatioFilter (self, ints, windowsize, damp):
    b = self._window(ints,windowsize,self._energy)
    B = self._windowStart(ints,self._energy)
    return b/(B+damp)

  @cython.wraparound(False)
  @cython.boundscheck(False)
  def energyRatio (self, traces=None, windowsize=earwindow, damp=eardamp):
    '''
    Return the energy ratio between a sliding window of a given size and the
    total energy to that point in the file.
    '''

    if (_isDescendant(traces, 'ndarray')):
      traces = traces.copy()
    else:
      traces = self.readTraces(traces)

    if (len(traces.shape) == 1):
      traces = self._energyRatioFilter(traces, windowsize, damp)
    else:
      for i in xrange(traces.shape[0]):
        traces[i][:] = self._energyRatioFilter(traces[i], windowsize, damp)[:]

    return traces
