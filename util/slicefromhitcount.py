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

import os
import warnings
import numpy as np
from optparse import OptionParser
from pygeo.fastio import getParams, readfast, writefast
from scipy.interpolate import griddata

# ------------------------------------------------------------------------
# Settings

AUTHORSHIP = 'Brendan Smithyman'
VERSION = '%prog v1.0\n'
DESCRIPTION = 'Uses ray hit count information to compute a 2D FAST model from the weighted average of cells in a 3D FAST model.'
USAGE = '%prog [options] input [output]'

addfilename = '.2D'
subsampleinterval = 10

# ------------------------------------------------------------------------

warnings.simplefilter('ignore', RuntimeWarning)

parser = OptionParser(	usage		= USAGE,
			version		= VERSION,
			description	= DESCRIPTION)

parser.add_option('-n', '--hitcount', action='store', dest='hitcount',
		help='file containing hitcount information [%default]')

parser.add_option('-d', '--direct2d', action='store_true', dest='direct2d',
		help='directly generate a 2D FAST model [%default]')

parser.add_option('-e', '--extrapolate', action='store_true', dest='extrap',
		help='use 2D spline fitting to extrapolate outside the area with ray coverage [%default]')

parser.set_defaults(	hitcount	= 'num.cell',
			direct2d	= False,
			extrap		= False)

(options, args) = parser.parse_args()

hcfile = options.hitcount

if (len(args) < 1):
  parser.error('Please specify an input model name.')

infile = args[0]

if (not os.path.isfile(infile)):
  parser.error('Input file "%s" does not exist!'%(infile,))

if (not os.path.isfile(hcfile)):
  parser.error('Hit count file "%s" does not exist!'%(hcfile,))

if (len(args) < 2):
  outfile = infile + addfilename
else:
  outfile = args[1]

try:
  params = getParams()
except IOError:
  parser.error('FAST header file "for.header" does not exist; are you executing in the FAST working directory?')

dims = params['dims']
pseudodims = (dims[0],1,dims[2])

try:
  model3d = readfast(infile, dims)
except:
  parser.error('Problem reading "%s"!'%(infile,))

try:
  hitcount = readfast(hcfile, dims)
except:
  parser.error('Problem reading "%s"!'%(hcfile,))

# ------------------------------------------------------------------------

# Calculate average of the 3D model along the y-axis
model2davg = model3d.mean(axis=1)

# Calculate 2D hitcount by summing along the y-axis
hitcount2dsum = hitcount.sum(axis=1)

# Make a copy of the 2D hitcount and reproject into 3D filled matrix
hitcount3dsum = hitcount2dsum.copy().reshape(pseudodims)

# Take the weighted average of the cells in each x,z location along the
# y-axis, normalized by the total number of rays in that x,z location
avg2d = np.nan_to_num(((1. * hitcount * model3d) / hitcount3dsum).sum(axis=1))

# Create mask (value of 1 where there are no rays)
selector = hitcount2dsum > 0
baseblank = 1 - selector

if (options.extrap):
  nx, nz = avg2d.shape
  X, Z = np.mgrid[0:nx, 0:nz]

  x0 = X[selector]
  z0 = Z[selector]
  pos0 = np.array([x0.ravel(),z0.ravel()]).T
  val = avg2d[selector]
  
  outim = griddata(pos0, val, P, method='nearest').reshape((nx,nz))

  # Locate missing points and extrapolate using bisplev
  filllocs = np.argwhere(baseblank)
  ravelfilllocs = np.argwhere(baseblank.ravel())
  X = filllocs[:,0]
  Y = filllocs[:,1]
  Z = interpolate.bisplev(X, Y, terms)
  model2d[ravelfilllocs] = Z

else:
  # Create an empty (2D) model in which to store the result
  # Copy the 2D average model into the masked region
  model2d = avg2d.copy() + model2davg.copy() * baseblank

if (options.direct2d):
  writefast(outfile, model2d.reshape(pseudodims))
else:
  newmodel3d = model3d.copy()
  newmodel3d[:] = model2d.reshape(pseudodims)
  writefast(outfile, newmodel3d)
