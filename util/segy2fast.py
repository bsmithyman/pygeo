#!/usr/bin/env python

import numpy as np
import sys
from optparse import OptionParser
from pygeo.segyread import SEGYFile
from pygeo.coord import reduceToLocal

# ------------------------------------------------------------------------
# Settings

AUTHORSHIP = 'Brendan Smithyman'
VERSION = '%prog v1.1\n'
DESCRIPTION = 'Exports a series of FAST pick datafiles based on SEG-Y headers.'
USAGE = '%prog [options] segy_file'


format_string = '%10.3f%10.3f%10.3f%10.3f%10.3f%3d\n'
unit = 1e3 # kilometres
zantithesis = -1
tfac = 1e3 # ms / sec
pickkey = 'delrt'
shotout = 'shotout.dat'
error = 30.

# ------------------------------------------------------------------------

parser = OptionParser(	usage		= USAGE,
			version		= VERSION,
			description	= DESCRIPTION)

parser.add_option('-b', '--basis', action='store', dest='basis',
		help='point to use as zero coordinate [%default]')

parser.add_option('-a', '--angle', action='store', dest='angle',
		help='angle in degrees for coordinate rotation [%default]')

parser.set_defaults(	basis = '0.,0.,0.',
			angle = '0.')

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


# Open SEG-Y file and get first trace header
sys.stdout.write('Reading "%s"...\n'%(infile,))
sys.stdout.flush()
sf = SEGYFile(infile, endian='Big')
trh0 = sf.trhead[0]

sys.stdout.write('Calculating scale factors...\n')
sys.stdout.flush()

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

sys.stdout.write('Writing output files...\n')
sys.stdout.flush()

# Create some empty lists to hold upcoming values
shotlocs = [[],[],[]]
shotactive = []

# Create bound thresholds (which will be updated)
bounds = [1e10,-1e10,1e10,-1e10,1e10,-1e10]

# Loop over each shot gather
for i in xrange(ngathers):

  # Create a FAST output file for this gather (using 4-digit filenames)
  outfile = 'fd%04d.ascii'%(i+1,)
  with open(outfile, 'w') as fp:
    sys.stdout.write('%s <-- SHOTID %d\n'%(outfile, shotnums[i]))
    sys.stdout.flush()

    # Get the trace header for the first trace in this shot gather
    trhl0 = sf.trhead[sf.ensembles[shotnums[i]]]
    sx = trhl0['sx'] * scalco
    sy = trhl0['sy'] * scalco
    sz = trhl0['selev'] * scalel * zantithesis

    (nsx, nsy, nsz) = reduceToLocal(np.array([sx,sy,sz],ndmin=2), angle, basis)[0]

    # Append information about this shot to the running tally of all shot
    # locations; this is used to construct f.in
    shotlocs[0].append(nsx)
    shotlocs[1].append(nsy)
    shotlocs[2].append(nsz)

    fp.write(format_string % (nsx, nsy, nsz, 0., 0., -1))

    tr0 = sf.ensembles[shotnums[i]]

    if (i == ngathers - 1):
      tr1 = sf.ntr - 1
    else:
      tr1 = sf.ensembles[shotnums[i+1]]

    shotactive.append(0)

    for j in xrange(tr0, tr1):
      trhl = sf.trhead[j]
      rx = trhl['gx'] * scalco
      ry = trhl['gy'] * scalco
      rz = trhl['gelev'] * scalel * zantithesis

      (nrx, nry, nrz) = reduceToLocal(np.array([rx,ry,rz],ndmin=2), angle, basis)[0]

      if (nrx < bounds[0]):
        bounds[0] = nrx
      if (nrx > bounds[1]):
        bounds[1] = nrx
      if (nry < bounds[2]):
        bounds[2] = nry
      if (nry > bounds[3]):
        bounds[3] = nry
      if (nrz < bounds[4]):
        bounds[4] = nrz
      if (nrz > bounds[5]):
        bounds[5] = nrz
      
      stime = trhl[pickkey]
      if ((stime != 0) and (stime != 65535)):
        fp.write(format_string % (nrx, nry, nrz, stime/tfac, error/tfac, 1))
        shotactive[-1] += 1

itrace = []
for i in xrange(ngathers):
  if (shotactive[i] > 0):
    itrace.append(i + 1)

with open(shotout, 'w') as fp:
  fp.write('      isource=')
  fp.write(', '.join(['%d'%(item != 0) for item in shotactive]))
  fp.write(',\n      xsource=')
  fp.write(', '.join(['%8.3f'%item for item in shotlocs[0]]))
  fp.write(',\n      ysource=')
  fp.write(', '.join(['%8.3f'%item for item in shotlocs[1]]))
  fp.write(',\n      zsource=')
  fp.write(', '.join(['%8.3f'%item for item in shotlocs[2]]))
  fp.write(',\n      itrace=')
  fp.write(', '.join(['%d'%item for item in itrace]))

sys.stdout.write('\nBounds:\n\t%f < x < %f\n\t%f < y < %f\n\t%f < z < %f\n' % tuple(bounds))
sys.stdout.flush()

