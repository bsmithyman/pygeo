
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

def _blockdims (dims, tilesize):
  return [item if (item)%2 else item+1 for i,item in enumerate(list(np.array(dims) / np.array(tilesize) + 1))]

def chequerboard (dims, tilesize=1, tilemin=0, tilemax=1):
  '''
  Creates a chequerboard ndarray with values of either tilemin or tilemax 
  (i.e. [0,1] by default). By default the tilesize one cell.
  '''

  if (not np.iterable(tilesize)):
    tilesize = tuple([tilesize for item in xrange(len(dims))])

  if (len(dims) != len(tilesize)):
    raise Exception()

  bd = _blockdims(dims, tilesize)

  a = 0
  b = 1
  for i in xrange(len(dims)):
    ind = -(i+1)
    vec = [tilesize[ind]*[a] if bit%2 else tilesize[ind]*[b] for bit in range(bd[ind])]
    a = np.array(vec)
    b = 1 - a

  vec = a
  vec = vec.reshape(tuple(np.array(bd) * np.array(tilesize)))
  vec = vec[[slice(0,endpoint) for endpoint in dims]]

  vec = vec * (tilemax - tilemin) + tilemin

  return vec

