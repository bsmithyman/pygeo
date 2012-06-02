#!/usr/bin/env python
from pygeo.autopick import energyRatio, energyRatioCL
from pygeo.segyread import SEGYFile
import pylab
import sys

infile = '/mnt/smaug/seisdata/line5-postpick.sgy'#'/spin/seisdata/line15-postpick.sgy'

def printnow (text):
  sys.stdout.write(text + '\n')
  sys.stdout.flush()

printnow('Reading SEGY data from %s'%(infile,))
sf = SEGYFile(infile, endian='Big')
a = sf.sNormalize(sf[:7550]) # a bit over 60 MB of data

printnow('Computing STA/LTA:')

printnow('\twith OpenMP...')
b = energyRatio(a)

printnow('\twith OpenCL...')
c = energyRatioCL(a)

printnow('\tDone.')


printnow('Plotting...')

pylab.figure()

pylab.subplot(3,1,1)
pylab.imshow(a.T)
pylab.title('Original Waveforms')

pylab.subplot(3,1,2)
pylab.imshow(b.T)
pylab.title('STA/LTA with OpenMP')

pylab.subplot(3,1,3)
pylab.imshow(c.T)
pylab.title('STA/LTA with OpenCL')

printnow('\tDone.')
pylab.show()
