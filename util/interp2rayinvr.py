#!/usr/bin/env python

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
from scipy.interpolate import griddata
from pygeo.fastio import readfast

infile = 'vel10.2D'
outfile = 'line15.rayinvr'

y = 0.
olddims = (151,65,61)
newdims = (301,121)
oldori = (-1.,-1.1,-2000.)
newori = (-1.,-2000.)
sl = 20
dxold = 100
dxnew = 50
forstr = '%f\t%f\t%f\t%f\n'

oxgrid, ozgrid = np.mgrid[oldori[0]:oldori[0]+olddims[0]*dxold:dxold,
                          oldori[2]:oldori[2]+olddims[2]*dxold:dxold]

nxgrid, nzgrid = np.mgrid[newori[0]:newori[0]+newdims[0]*dxnew:dxnew,
                          newori[1]:newori[1]+newdims[1]*dxnew:dxnew]

fmodel = readfast(infile, olddims)
oldsl = fmodel[:,sl,:]

gridinput = np.column_stack((oxgrid.ravel(), ozgrid.ravel()))
gridoutput = np.column_stack((nxgrid.ravel(), nzgrid.ravel()))

newsl = griddata(gridinput, oldsl.ravel(), gridoutput, method='linear')
newsl.shape = newdims

saveoutput = np.column_stack((nxgrid.T.ravel(), nzgrid.T.ravel(),newsl.T.ravel()))

with open(outfile, 'w') as fp:
  for i in range(saveoutput.shape[0]):
    fp.write(forstr%(saveoutput[i,0], y, -saveoutput[i,1], saveoutput[i,2]))
