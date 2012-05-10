#!/usr/bin/env python

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

if (len(args < 2)):
  parser.print_help()
  exit(0)

infile = args[0]
outfile = args[1]

if (options.dims is None):
  dims = getParams()['dims'] 
else:
  dims = [int(item) for item in options.dims.strip().split(',')]

velmodel = readfast(infile, dims).copy(order='F')
with open(modelout, 'w') as fp:
  fp.write(velmodel)
