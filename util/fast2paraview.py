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
from pygeo.fastio import readfast, getParams
from optparse import OptionParser

AUTHORSHIP = 'Brendan Smithyman, May 2012'
VERSION = '%prog v1.0\n'
DESCRIPTION = '''Convert FAST model to Paraview RAW file.'''
USAGE = '%prog [options] infile outfile'

parser = OptionParser(	usage		= USAGE,
			version		= VERSION,
			description	= DESCRIPTION)

parser.add_option('-d', '--dims', action='store', dest='dims',
		help='dimensions of FAST model')

parser.set_defaults(	dims = None)

(options, args) = parser.parse_args()

if (len(args) < 2):
  parser.print_help()
  exit(0)

infile = args[0]
outfile = args[1]

if (options.dims is None):
  dims = getParams()['dims'] 
else:
  dims = [int(item) for item in options.dims.strip().split(',')]

velmodel = readfast(infile, dims).copy(order='F')
with open(outfile, 'w') as fp:
  fp.write(velmodel)
