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

import sys
import numpy as np

if (len(sys.argv) != 3):
  print('Run as:\n\trinigengroups.py projnm nsg')
  exit()

projnm, nsg = sys.argv[1:]
nsg = int(nsg)

header1 = '<isg ><isg ><isg ><isg ><isg ><isg ><isg ><isg >\n'
header2 = '<dsd   ><dsd   ><dsd   ><dsd   ><dsd   ><dsd   >\n'

data1part = '% 6d'
data1 = data1part*8 + '\n'

data2part = '% 8.3f'
data2 = data2part*6 + '\n'

outlines = []

for i in xrange(nsg/8):
  chunk = data1 % tuple(range(i*8 + 1, i*8 + 9))
  sys.stdout.write(header1)
  sys.stdout.write(chunk)
  outlines.append(header1)
  outlines.append(chunk)

sys.stdout.write(header1)
outlines.append(header1)
outlines.append('')

nsgremainder = np.mod(nsg, 8)
nsgbase = nsg - nsgremainder
for i in xrange(nsgremainder):
  chunk = data1part % (nsgbase + i + 1)
  sys.stdout.write(chunk)
  outlines[-1] += chunk
outlines[-1] += '\n'
sys.stdout.write('\n')

for i in xrange(nsg/6):
  chunk = data2 % tuple([1.]*6)
  sys.stdout.write(header2)
  sys.stdout.write(chunk)
  outlines.append(header2)
  outlines.append(chunk)

sys.stdout.write(header2)
outlines.append(header2)
outlines.append('')

nsgremainder = np.mod(nsg, 6)
nsgbase = nsg - nsgremainder
for i in xrange(nsgremainder):
  chunk = data2part % (1.,)
  sys.stdout.write(chunk)
  outlines[-1] += chunk
outlines[-1] += '\n'
sys.stdout.write('\n')

f = open(projnm + '.partrini', 'w')
f.writelines(outlines)
f.close()
