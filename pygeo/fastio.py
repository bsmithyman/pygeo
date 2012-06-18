
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
import mmap
import struct

rlpad = 4

if (rlpad == 8):
  rlpak = 'Q'
elif (rlpad == 4):
  rlpak = 'I'
else:
  raise FIOUnrecognizedRecordLength

# ----------------------------------------------------------------------
def readfast (filename, dims):
  '''
  Read FAST model file to Python ndarray.
  '''

  (nx, ny, nz) = dims

  fp = open(filename, 'rb')
  map = mmap.mmap(fp.fileno(), 0, flags=mmap.MAP_PRIVATE)

  vol = np.zeros(dims, dtype=np.int16)
  for k in range(nz):
    for j in range(ny):
      bytestart = k*ny*(2*rlpad+2*nx) + j*(2*rlpad+2*nx)
      entries = struct.unpack('%dH'%(nx,), map[bytestart+rlpad:bytestart+rlpad+nx*2])
      vol[:,j,k] = np.array(entries)

  map.close()
  fp.close()

  return vol

# ----------------------------------------------------------------------
def writefast (filename, model):
  '''
  Write FAST model file from Python ndarray.
  '''

  (nx, ny, nz) = model.shape
  rmark = struct.pack(rlpak, np.int64(nx*2))

  fp = open(filename, 'wb')

  for k in range(nz):
    for j in range(ny):
      outstring = struct.pack('%dH'%(nx,), *list(model[:,j,k]))
      fp.write(rmark)
      fp.write(outstring)
      fp.write(rmark)

  fp.close()

# ----------------------------------------------------------------------
def readbath (filename, dims):
  '''
  Read FAST bathymetry file to Python ndarray.
  '''
  
  (nx, ny, nz) = dims

  fp = open(filename, 'rb')
  map = mmap.mmap(fp.fileno(), 0, flags=mmap.MAP_PRIVATE)

  bathy = np.zeros(dims[:2], dtype=np.float32)
  for j in range(ny):
    bytestart = j*(4*nx+2*rlpad)
    bathy[:,j] = np.array(struct.unpack('%df'%(nx,), map[bytestart+rlpad:bytestart+rlpad+4*nx]))

  map.close()
  fp.close()

  return bathy

# ----------------------------------------------------------------------
def restorebath (modeltop, modelbottom, bathymetry, dims, bounds):
  '''
  Restore top of the model above the bathymetry level.
  '''

  (nx, ny, nz) = dims
  zspace = np.linspace(bounds[2][0],bounds[2][1],dims[2])
  
  result = modelbottom.copy()
  
  for i in range(dims[0]):
    for j in range(dims[1]):
      for k in range(dims[2]):
        if (zspace[k] <= bathymetry[i][j]):
          result[i,j,k] = modeltop[i,j,k]
      #result[i,j,:] = [v for v in modeltop[i,j,:] if 
  return result

# ----------------------------------------------------------------------
def readPicks (filename):

  with open(filename, 'r') as fp:
    lines = fp.readlines()

  block = np.array([[float(item) for item in line.strip().split()] for line in lines])

  sourcecoords = block[0,:3]
  reccoords = block[1:,:3]
  times = block[1:,3]
  errors = block[1:,4]

  return {'sc': sourcecoords, 'rc': reccoords, 't': times, 'e': errors}

# ----------------------------------------------------------------------
def getParams (filename = 'for.header'):
  with open(filename, 'r') as fp:
    res = fp.read()

  res = res.strip().split()
  paramdict = {}
  paramdict['x'] = [float(item) for item in res[0:2]]
  paramdict['y'] = [float(item) for item in res[2:4]]
  paramdict['z'] = [float(item) for item in res[4:6]]
  paramdict['dx'] = float(res[6])
  paramdict['dims'] = [int(item) for item in res[7:10]]

  return paramdict
