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
