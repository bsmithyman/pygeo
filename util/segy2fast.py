#!/usr/bin/env python

from pygeo.segyread import SEGYFile
import sys

# ------------------------------------------------------------------------

format_string = '%10.3f%10.3f%10.3f%10.3f%10.3f%3d\n'
unit = 1e3 # kilometres
zantithesis = -1
tfac = 1e3 # ms / sec
pickkey = 'delrt'
shotout = 'shotout.dat'

# ------------------------------------------------------------------------

infile = sys.argv[1]

sys.stdout.write('Reading "%s"...\n'%(infile,))
sys.stdout.flush()
sf = SEGYFile(infile)
trh0 = sf.trhead[0]

sys.stdout.write('Calculating scale factors...\n')
sys.stdout.flush()

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

sys.stdout.write('Calculating ensemble boundaries...\n')
sys.stdout.flush()
sf._calcEnsembles()

ngathers = len(sf.ensembles)
shotnums = sf.ensembles.keys()
shotnums.sort()

sys.stdout.write('Writing output files...\n')
sys.stdout.flush()

shotlocs = [[],[],[]]
shotactive = []

bounds = [1e10,-1e10,1e10,-1e10,1e10,-1e10]
for i in xrange(ngathers):
  outfile = 'fd%04d.ascii'%(i+1,)
  with open(outfile, 'w') as fp:
    sys.stdout.write('%s <-- SHOTID %d\n'%(outfile, shotnums[i]))
    sys.stdout.flush()
    trhl0 = sf.trhead[sf.ensembles[shotnums[i]]]
    sx = trhl0['sx'] * scalco
    sy = trhl0['sy'] * scalco
    sz = trhl0['selev'] * scalel * zantithesis
    shotlocs[0].append(sx)
    shotlocs[1].append(sy)
    shotlocs[2].append(sz)

    fp.write(format_string % (sx, sy, sz, 0., 0., -1))

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

      if (rx < bounds[0]):
        bounds[0] = rx
      if (rx > bounds[1]):
        bounds[1] = rx
      if (ry < bounds[2]):
        bounds[2] = ry
      if (ry > bounds[3]):
        bounds[3] = ry
      if (rz < bounds[4]):
        bounds[4] = rz
      if (rz > bounds[5]):
        bounds[5] = rz
      
      stime = trhl[pickkey]
      error = 0.
      if ((stime != 0) and (stime != 65535)):
        fp.write(format_string % (rx, ry, rz, stime/tfac, error/tfac, 1))
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

