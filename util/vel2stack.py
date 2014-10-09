#!/usr/bin/env python

# pygeo - a distribution of tools for managing geophysical data
# Copyright (C) 2013 Brendan Smithyman

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

import os
import numpy as np
from pygeo.segyread import SEGYFile
from pygeo.coord import reduceToLocal

# ------------------------------------------------------------------------
# Settings

AUTHORSHIP = 'Brendan Smithyman, May 2013'
USAGE = '%prog [options] stackfile velfile outfile'
VERSION = '%prog v1.0\n'
DESCRIPTION = '''
Convert 2D velocity model to CDP-located stacking-style image.
'''

# ------------------------------------------------------------------------
# Parameter parsing

from optparse import OptionParser

parser = OptionParser(	usage		= USAGE,
			version		= VERSION,
			description	= DESCRIPTION)

parser.add_option('-b', '--basis', action='store', dest='basis',
		help='point to use as zero coordinate [%default]')

parser.add_option('-a', '--angle', action='store', dest='angle',
		help='angle in degrees for coordinate rotation [%default]')

parser.add_option('-u', '--unit', action='store', dest='unit',
		help='spatial unit [%default]')

parser.add_option('-z', '--zdir', action='store', dest='zdir',
		help='coord. system z-scaling [%default]')

parser.add_option('-e', '--extent', action='store', dest='extent',
		help='extent of velocity model (x0,x1) [%default]')

parser.set_defaults(	basis	= '0.,0.,0.',
			angle	= '0.',
			unit	= '1e0',
			zdir	= '-1',
			extent = None,
)

(options, args) = parser.parse_args()

if (len(args) < 3):
  parser.error('Please specify an input stack and velocity model, and an output file!')

if (not options.extent):
  parser.error('Please specify the extent of the velocity model!')

x0, x1 =  [np.float(item) for item in options.extent.strip().split(',')]

# Convert rotation angle to radians
angle = np.float(options.angle)*np.pi/180.

# Convert basis to array
basis = np.array([np.float(item) for item in options.basis.strip().split(',')])

unit = np.float(options.unit)
zantithesis = np.float(options.zdir)

stackfile = args[0]
velfile = args[1]
outfile = args[2]

sfstack = SEGYFile(stackfile, endian='Big')
sfvel = SEGYFile(velfile, endian='Big')

# ------------------------------------------------------------------------
# Main Program

strh0 = sfstack.trhead[0]
sntr = sfstack.ntr
sns = sfstack.ns

vtrh0 = sfvel.trhead[0]
vntr = sfvel.ntr
vns = sfvel.ns

scalco = float(strh0['scalco'])
if (scalco < 0):
  scalco = -1./scalco

scalel = float(strh0['scalel'])
if (scalel < 0):
  scalel = -1./scalel

coordarr = np.zeros((sntr,3), dtype=np.float32)

coordarr[:,0] = np.array([float(trh['sx'])*scalco for trh in sfstack.trhead])
coordarr[:,1] = np.array([float(trh['sy'])*scalco for trh in sfstack.trhead])

newcoordarr = reduceToLocal(coordarr, angle, basis)

velarr = sfvel[:]
outarr = np.zeros((sntr,sns))

vdx = (x1 - x0) / (vntr - 1)

for i in xrange(sntr):
  newx = newcoordarr[i, 0]
  vi = np.round((newx - x0) / vdx)
  outarr[i,:] = velarr[vi,:]

sfstack.writeSU(outfile + '.su', outarr)
os.system('< %s segyhdrs | segywrite endian=0 tape=%s'%(outfile + '.su', outfile))
