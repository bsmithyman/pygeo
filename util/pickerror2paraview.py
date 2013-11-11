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

import glob
import numpy as np
import pygeo.fastio as fastio
import vtk

origfilter = '*.ascii'
calcfilter = '*.calcascii'
outfile = 'errors.vtp'
coord_scale = 1000.
vertical_scale = 1000.

origfiles = glob.glob(origfilter)
calcfiles = glob.glob(calcfilter)

origfiles.sort()
calcfiles.sort()

maxiter = max(len(calcfiles),len(origfiles))
miniter = min(len(calcfiles),len(origfiles))


if (miniter != maxiter):
  print('Working from a subset of all shot gathers: %d/%d'%(miniter,maxiter))
  print('At present this doesn\'t work!')
  numgathers = miniter
  exit(1)
    
else:
  numgathers = maxiter

midpoints = []
offsets = []
errors = []
relerrors = []

for i in xrange(numgathers):
  orig = fastio.readPicks(origfiles[i])
  calc = fastio.readPicks(calcfiles[i])

  midpoints.append((orig['rc'] + orig['sc'])/2.)
  offsets.append(orig['rc'] - orig['sc'])
  errors.append(calc['t'] - orig['t'])
  relerrors.append(errors[-1] / orig['e'])

midpoints = np.concatenate(midpoints)
offsets = np.concatenate(offsets)
errors = np.concatenate(errors)
relerrors = np.concatenate(relerrors)

# ------------------------------------------------------------------------

outcoords = midpoints.copy()
outcoords = outcoords * coord_scale
absoffsets = np.sqrt((offsets**2).sum(axis=1))
absoffsets = absoffsets * vertical_scale
outcoords[:,2] = -absoffsets

coordarray = outcoords.ravel()
valarray = errors.ravel()

points = vtk.vtkPoints()
darr = vtk.vtkDoubleArray()
darr.SetVoidArray(coordarray, len(coordarray), 1)
darr.SetNumberOfComponents(3)
points.SetData(darr)

vals = vtk.vtkDoubleArray()
vals.SetVoidArray(valarray, len(valarray), 1)
vals.SetNumberOfComponents(1)
vals.SetName('error')

polydata = vtk.vtkPolyData()
polydata.SetPoints(points)
polydata.GetPointData().SetScalars(vals)

writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(outfile)
writer.SetInput(polydata)
writer.Write()
