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
import matplotlib.pylab as pl
from pygeo.fullpy import readini
from pygeo.segyread import SEGYFile

# ------------------------------------------------------------------------
# Settings

AUTHORSHIP = 'Brendan Smithyman, October 2012'
VERSION = '%prog v1.0\n'
DESCRIPTION = '''
Generates model plots with useful information overlaid; outputs guideline
information for troubleshooting.
'''

# ------------------------------------------------------------------------
# Parameter parsing

from optparse import OptionParser

usage = '%prog [options] [projectname]'

parser = OptionParser(usage       = usage,
                      version     = VERSION,
                      description = DESCRIPTION)

parser.add_option('-v', '--verbose', action='store_true', dest='verbose',
		help='display additional information')

parser.set_defaults(
			verbose		= False,
)

(options, args) = parser.parse_args()

if (len(args) == 0):
  parser.error('Please specify a filename!')

gopts = {}

gopts['filename'] = args[0]

try:
  globresult = glob.glob('*.ini')
  if (len(globresult) != 0):
    for gri in globresult:
      thisprojnm = gri[:-4]
      if (gopts['filename'].find(thisprojnm) == 0):
        gopts['projnm'] = thisprojnm
        gopts['inifile'] = gri
        if (options.verbose):
          print('Found matching project file: %s'%(gri,))
        break
  else:
    parser.error('This program requires a project *.ini file to function!')
except:
  parser.error('Something happened!')

ini = readini(gopts['inifile'])

extent = [ini['xorig'],
          ini['xorig'] + ini['dx']*(ini['nx']-1),
        -(ini['zorig'] + ini['dz']*(ini['nz']-1)),
         -ini['zorig']]

srcs = ini['srcs']

sf = SEGYFile(gopts['filename'], endian='Big', verbose=options.verbose)
model = sf[:].T

fig = pl.figure()
ax = fig.add_subplot(1,1,1, aspect=1.0)
pl.imshow(model, extent=extent)
pl.plot(srcs[:,0], -srcs[:,2], 'w.')
ax.axis(extent)

pl.show()
