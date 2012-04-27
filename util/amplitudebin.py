#!/usr/bin/env python
import numpy as np
import cPickle as pickle
import sys

from pygeo.segyread import *

infile = sys.argv[1]
outfile = sys.argv[2]
scalco = float(sys.argv[3])

def offsets (sf):
  result = np.empty((len(sf),))
  for i in range(len(sf)):
    h = sf.trhead[i]
    result[i] = np.sqrt((h['gx']-h['sx'])**2 + (h['gy']-h['sy'])**2 + (h['gelev'] - h['selev'])**2)

  return result*scalco

def amplitudes (sf):
  result = np.empty((len(sf),))
  for i in range(len(sf)):
    result[i] = np.sqrt((sf.readTraces(i+1)**2).sum())

  return result

print('Reading SEGY file [%s]...'%infile)
sf = SEGYFile(infile)

print('Computing offsets using scale of [%f]...'%scalco)
offs = offsets(sf)

print('Computing amplitudes...')
amps = amplitudes(sf)

print('Saving results [%s]...'%outfile)
with open(outfile, 'wb') as fp:
  cp = pickle.Pickler(fp)
  cp.dump(offs)
  cp.dump(amps)
  del cp

