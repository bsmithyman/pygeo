#!/usr/bin/env python

import numpy as np
from pygeo.segyread import SEGYFile
import sys

datum = (0.,0.,0.)
outfiles = ['shotgeom.txt', 'recgeom.txt']

forma = '%8d % 10.5E % 10.5E % 10.5E %7.3f\n'
zantithesis=1

infile = sys.argv[1]

if (len(sys.argv) == 5):
  datum = (float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]))

sf = SEGYFile(infile, endian='Big')
sf._calcEnsembles()
trh0 = sf.trhead[0]

if (trh0['scalco'] < 0):
  scalco = 1./abs(trh0['scalco'])
else:
  scalco = trh0['scalco']

if (trh0['scalel'] < 0):
  scalel = 1./abs(trh0['scalel'])
else:
  scalel = trh0['scalel']

rectrh = np.array([(float(trh['gx']),float(trh['gy']),float(trh['gelev'])) for trh in sf.trhead[:755]], dtype=np.float32)*scalco

reclines = []
for i in xrange(len(rectrh)):
  reclines.append(forma%(i+1,rectrh[i,0] - datum[0],rectrh[i,1] - datum[1],rectrh[i,2] - datum[2],1.))

ngathers = len(sf.ensembles)
ordering = np.argsort(sf.ensembles.values())
shotnums = np.array(sf.ensembles.keys())[ordering]

shotlines = []
for i in xrange(ngathers):
  trhl0 = sf.trhead[sf.ensembles[shotnums[i]]]
  sx = trhl0['sx'] * scalco
  sy = trhl0['sy'] * scalco
  sz = trhl0['selev'] * scalel * zantithesis
  shotlines.append(forma%(i+1,sx - datum[0],sy - datum[1],sz - datum[2],1.))

with open(outfiles[0], 'w') as fp:
  fp.writelines(shotlines)

with open(outfiles[1], 'w') as fp:
  fp.writelines(reclines)
