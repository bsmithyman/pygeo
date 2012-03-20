from scipy import zeros, array, int16, int64, float32, linspace
from mmap import *
from struct import *

rlpad = 4

if (rlpad == 8):
  rlpak = 'Q'
elif (rlpad == 4):
  rlpak = 'L'
else:
  raise FIOUnrecognizedRecordLength

# ----------------------------------------------------------------------
def readfast (filename, dims):
  '''
  Read FAST model file to Python ndarray.
  '''

  (nx, ny, nz) = dims

  fp = open(filename, 'rb')
  map = mmap(fp.fileno(), 0, flags=MAP_PRIVATE)

  vol = zeros(dims, dtype=int16)
  for k in range(nz):
    for j in range(ny):
      bytestart = k*ny*(2*rlpad+2*nx) + j*(2*rlpad+2*nx)
      entries = unpack('%dH'%(nx,), map[bytestart+rlpad:bytestart+rlpad+nx*2])
      vol[:,j,k] = array(entries)

  map.close()
  fp.close()

  return vol

# ----------------------------------------------------------------------
def writefast (filename, model):
  '''
  Write FAST model file from Python ndarray.
  '''

  (nx, ny, nz) = model.shape
  rmark = pack(rlpak, int64(nx*2))

  fp = open(filename, 'wb')

  for k in range(nz):
    for j in range(ny):
      outstring = pack('%dH'%(nx,), *list(model[:,j,k]))
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
  map = mmap(fp.fileno(), 0, flags=MAP_PRIVATE)

  bathy = zeros(dims[:2], dtype=float32)
  for j in range(ny):
    bytestart = j*(4*nx+2*rlpad)
    bathy[:,j] = array(unpack('%df'%(nx,), map[bytestart+rlpad:bytestart+rlpad+4*nx]))

  map.close()
  fp.close()

  return bathy

# ----------------------------------------------------------------------
def restorebath (modeltop, modelbottom, bathymetry, dims, bounds):
  '''
  Restore top of the model above the bathymetry level.
  '''

  (nx, ny, nz) = dims
  zspace = linspace(bounds[2][0],bounds[2][1],dims[2])
  
  result = modelbottom.copy()
  
  for i in range(dims[0]):
    for j in range(dims[1]):
      for k in range(dims[2]):
        if (zspace[k] <= bathymetry[i][j]):
          result[i,j,k] = modeltop[i,j,k]
      #result[i,j,:] = [v for v in modeltop[i,j,:] if 
  return result
