#!/usr/bin/env python2

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

