#!/usr/bin/env python

import numpy as np
from pygeo.segyread import SEGYFile
from pygeo.coord import reduceToLocal
import sys
from optparse import OptionParser

# ------------------------------------------------------------------------
# Settings

AUTHORSHIP = 'Brendan Smithyman'
VERSION = '%prog v1.1\n'
DESCRIPTION = 'Exports a series of FULLWV geometry files based on SEG-Y headers.'
USAGE = '%prog [options] segy_file'

output_format = '%8d % 10.5E % 10.5E % 10.5E %7.3f\n'

# ------------------------------------------------------------------------

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

parser.add_option('-s', '--shotout', action='store', dest='shotout',
		help='filename for shot geometry [%default]')

parser.add_option('-r', '--recout', action='store', dest='recout',
		help='filename for receiver geometry [%default]')

parser.set_defaults(	basis	= '0.,0.,0.',
			angle	= '0.',
			unit	= '1e0',
			zdir	= '-1',
			shotout	= 'shotgeom.txt',
			recout	= 'recgeom.txt')

(options, args) = parser.parse_args()

if (len(args) < 1):
  parser.error('Please specify a SEG-Y file!')
  exit(1)

# Get input filename
infile = args[0]

# Convert rotation angle to radians
angle = np.float(options.angle)*np.pi/180.

# Convert basis to array
basis = np.array([np.float(item) for item in options.basis.strip().split(',')])

unit = np.float(options.unit)
zantithesis = np.float(options.zdir)
shotout = options.shotout
recout = options.recout

# Open SEG-Y file and get first trace header
sys.stdout.write('Reading "%s"...\n'%(infile,))
sys.stdout.flush()
sf = SEGYFile(infile, endian='Big')
trh0 = sf.trhead[0]

# Determine coordinate and elevation scale factors from first trace header
# (assume same for all traces)
if (trh0['scalco'] < 0):
  scalco = 1./abs(trh0['scalco'])
else:
  scalco = trh0['scalco']

scalco = scalco / unit

if (trh0['scalel'] < 0):
  scalel = 1./abs(trh0['scalel'])
else:
  scalel = trh0['scalel']

scalel = scalel / unit

# Use SEGYFile internal to calculate shot-gather boundaries
sys.stdout.write('Calculating ensemble boundaries...\n')
sys.stdout.flush()
sf._calcEnsembles()

# Find the number of ensembles, and order them by occurrence in the SEG-Y file
ngathers = len(sf.ensembles)
ordering = np.argsort(sf.ensembles.values())
shotnums = np.array(sf.ensembles.keys())[ordering]

# Assume that SEG-Y file is set up for position-consistent receivers
nrecs = sf.ensembles[shotnums[1]]

rectrh = np.array([(float(trh['gx'])*scalco,float(trh['gy'])*scalco,zantithesis*float(trh['gelev'])*scalel) for trh in sf.trhead[:nrecs]], dtype=np.float32)
newrectrh = reduceToLocal(rectrh, angle, basis)

reclines = []
for i in xrange(len(rectrh)):
  reclines.append(output_format%(i+1,newrectrh[i,0],newrectrh[i,1],newrectrh[i,2],1.))

shottrh = np.array([(float(trh['sx'])*scalco,float(trh['sy'])*scalco,zantithesis*float(trh['selev'])*scalel) for trh in [sf.trhead[sf.ensembles[sn]] for sn in shotnums]])
newshottrh = reduceToLocal(shottrh, angle, basis)

shotlines = []
for i in xrange(ngathers):
  shotlines.append(output_format%(i+1,newshottrh[i,0],newshottrh[i,1],newshottrh[i,2],1.))

with open(shotout, 'w') as fp:
  fp.writelines(shotlines)

with open(recout, 'w') as fp:
  fp.writelines(reclines)
