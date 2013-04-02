#!/usr/bin/env python2

# pygeo - a distribution of tools for managing geophysical data
# Copyright (C) 2011-2013 Brendan Smithyman

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

# ------------------------------------------------------------------------

import sys
import warnings
import numpy as np
from scipy.interpolate import griddata
from pygeo.segyread import SEGYFile
from pygeo.coord import reduceToLocal
from pygeo.fastio import writefast

# ------------------------------------------------------------------------
# Settings

AUTHORSHIP = 'Brendan Smithyman, April 2013'
USAGE = '%prog [options] infile [outfile]'
VERSION = '%prog v1.0\n'
DESCRIPTION = '''
Project 3D seismic point information (e.g., images) to 2D plane.
'''

# ------------------------------------------------------------------------
# Useful Functions

def printnow (text):
  sys.stdout.write(text + '\n')
  sys.stdout.flush()

# ------------------------------------------------------------------------
# Parameter parsing

from optparse import OptionParser

parser = OptionParser(	usage		= USAGE,
			version		= VERSION,
			description	= DESCRIPTION)

parser.add_option('-a', '--angle', action='store', dest='angle',
		help='angle in degrees for coordinate rotation [%default]')

parser.add_option('-b', '--basis', action='store', dest='basis',
		help='point to use as zero coordinate [%default]')

parser.add_option('-c', '--cellsize', action='store', dest='cellsize',
		help='cell size for output model [%default]')

parser.add_option('-g', '--geom', action='store', dest='geom',
		help='file for geometry information [%default]')

parser.add_option('-o', '--outbox', action='store', dest='outbox',
		help='local geometry extent used for output [%default]')

parser.add_option('-u', '--unit', action='store', dest='unit',
		help='spatial unit [%default]')

parser.set_defaults(	angle 		= '0.',
			basis		= '0.,0.,0.',
			cellsize	= 'Auto',
			geom		= None,
			outbox 		= 'Auto',
			unit		= '1e0',
)

(options, args) = parser.parse_args()

if (len(args) < 1):
  parser.error('Please specify an input file!')

# Get input filename
infile = args[0]

# Generate output filename
if (len(args) >= 2):
  outfile = args[1]
else:
  outfile = infile + '.proj'

# Convert rotation angle to radians
angle = np.float(options.angle)*np.pi/180.

# Convert basis to array
basis = np.array([np.float(item) for item in options.basis.strip().split(',')])

# Set units
unit = np.float(options.unit)

printnow('Reading input files')

# Get data and geometry input files
sfdata = SEGYFile(infile)
if (not options.geom is None):
  geomfile = options.geom
  sfgeom = SEGYFile(geomfile)
else:
  geomfile = infile
  sfgeom = sfdata

printnow('Computing geometry')

# Grab headers
trh0 = sfgeom.trhead[0]
ntr = sfgeom.ntr
ns = sfgeom.ns

# Check to see if geometry and data files are the same shape
if (ntr != sfdata.ntr or ns != sfdata.ns):
  parser.error('Geometry file and input file must have the same number of samples and traces!')

# Compute cell size
dz = trh0['dt']/1000.
if (options.cellsize == 'Auto'):
  cellsize = dz
else:
  cellsize = int(options.cellsize)

# Compute coordinate scaling factors
scalco = float(trh0['scalco'])
if (scalco < 0):
  scalco = -1./scalco
scalco = scalco / unit

scalel = float(trh0['scalel'])
if (scalel < 0):
  scalel = -1./scalel
scalel = scalel / unit

# Get the x,y,z coordinates for the top of each trace
traceX = np.array([float(trh['sx'])*scalco for trh in sfgeom.trhead])
traceY = np.array([float(trh['sy'])*scalco for trh in sfgeom.trhead])
traceZ = np.array([float(trh['selev'])*scalel for trh in sfgeom.trhead])

# Some working arrays
theones = np.ones((ntr,ns))
theslope = -np.arange(ns) * dz
coordarray = np.empty((ntr,ns,3), dtype=np.float64)

# Create an array of coordinates of shape (ntr, ns, 3)
coordarray[:,:,0] = (traceX * theones.T).T
coordarray[:,:,1] = (traceY * theones.T).T
coordarray[:,:,2] = (traceZ * theones.T).T + theslope * theones

# Reproject to the new coordinate system by rotating around the z-axis
newcoordarray = reduceToLocal(coordarray.reshape((ntr*ns,3)), angle, basis).reshape((ntr,ns,3))

newxcoords = newcoordarray[:,:,0]
newzcoords = newcoordarray[:,:,2]

# Set up a bounding box for the original extent of the file and possibly a different
# bounding box for the output
interbox = [newxcoords.min(), newxcoords.max(), newzcoords.min(), newzcoords.max()]
if (options.outbox == 'Auto'):
  outbox = interbox
else:
  outbox = np.array([np.float(item) for item in options.outbox.strip().split(',')])

# Determine the model size for the new output
nx = round((outbox[1] - outbox[0]) / cellsize) + 1
nz = round((outbox[3] - outbox[2]) / cellsize) + 1

printnow('Resampling image')

# Get image data
im = sfdata[:]

# Set up arrays to hold reprojected image
interim = np.zeros((nx,nz), dtype=np.float32)
interimcount = np.zeros((nx,nz))

# Reproject the image to the new geometry, computing the new cells by simple weighted
# nearest neighbour (doesn't work for interpolation, just for downsampling)
posarr = newcoordarray.reshape((ntr*ns,3))
valarr = im.reshape((ntr*ns))
for i in xrange(len(posarr)):
  # For each position in the input find the intersecting cell and
  # add to the accumulator
  x, y, z = posarr[i]
  ix = round((x - outbox[0])/cellsize)
  iz = round((outbox[3] - z)/cellsize)
  # Ignore things that are outside the bounds
  if (ix >= 0 and ix < nx and iz >= 0 and iz < nz):
    interimcount[ix,iz] += 1
    interim[ix,iz] += valarr[i]

# Apply the weighting, but ignore division by zero
with warnings.catch_warnings():
  warnings.simplefilter("ignore")
  interim /= interimcount

# Generate interior model without nans
interim = np.nan_to_num(interim)

# Create selection array
selector = interimcount != 0

# Generate output grids
X, Z = np.mgrid[0:nx, 0:nz]
X = X*cellsize + outbox[0]
Z = outbox[3] - Z*cellsize
P = np.array([X.ravel(),Z.ravel()]).T

#Generate initial grids
x0 = X[selector]
z0 = Z[selector]
pos0 = np.array([x0.ravel(),z0.ravel()]).T
val = interim[selector]

printnow('Extrapolating')

# Extrapolate using nearest neighbour
outim = griddata(pos0, val, P, method='nearest').reshape((nx,nz))

printnow('Writing output')

# Write FAST model
outdims = (nx, 1, nz)
writefast(outfile, outim.reshape(outdims)
