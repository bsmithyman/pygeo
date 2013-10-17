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
import cPickle as pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm

al = [0,14400,1E-8,1E+4]
lal = [al[0],al[1],np.log(al[2]),np.log(al[3])]
binbounds = [25, 1000, 13500]
threshcnt = 10
filelist = sys.argv[1:3]#['orig.out', 'filtered.out']
eps = 1E-10
polyorder = 1

def bindata (offs, amps):
  [nbins, x0, x1] = binbounds
  valid = np.argwhere((amps > al[2])*(amps < al[3]))
  coffs = offs[valid]
  camps = amps[valid]
  bnds = np.linspace(x0, x1, nbins+1)
  accum = np.zeros((nbins,), dtype=np.float64)
  
  for i in xrange(nbins):
    indexarr = np.argwhere((coffs >= bnds[i])*(coffs < bnds[i+1]))
    accum[i] = camps[indexarr].mean()

  midpts = (bnds[1:]+bnds[:-1])/2

  return midpts, accum

reslist = []
for file in filelist:

  with open(file, 'rb') as fp:
    cp = pickle.Unpickler(fp)
    offs = cp.load()
    amps = cp.load()
    del cp

  lamps = np.log(amps)
  reslist.append([offs, amps, lamps])

ox, oy = bindata(reslist[0][0], reslist[0][1])
loy = np.log(oy)
pfo = np.polyfit(ox, loy, polyorder)
ploy = np.polyval(pfo, ox)

fx, fy = bindata(reslist[1][0], reslist[1][1])
lfy = np.log(fy)
pff = np.polyfit(fx, lfy, polyorder)
plfy = np.polyval(pff, fx)

diffy = loy - lfy
pdiffterm = np.polyfit(ox, diffy, polyorder)
pdiffy = np.polyval(pdiffterm, ox)


plt.figure()
plt.clf()
plt.hsv()

ax1 = plt.axes([0.125,0.1,0.775,0.15])
plt.plot(ox, diffy, 'kx')
plt.plot(ox, pdiffy, 'k-')
aln = plt.axis()
plt.axis([al[0],al[1],aln[2],aln[3]])
plt.annotate('$A_1 = %e$\n$A_0 = %e$'%tuple(-pdiffterm), [0.9*(al[1]-al[0])+aln[0], 0.4*(aln[3]-aln[2])+aln[2]])
plt.ylabel('$\Delta \ln (A)$')
plt.xlabel('Offset (m)')

ax2 = plt.axes([0.125,0.3,0.775,0.6], sharex=ax1)
plt.hexbin(reslist[0][0], reslist[0][2], mincnt=threshcnt)
plt.plot(ox, loy, 'kx')
plt.plot(ox, ploy, 'k-')
plt.annotate('Field', [0.95*al[1], loy.mean()])
plt.hexbin(reslist[1][0], reslist[1][2], mincnt=threshcnt)
plt.plot(fx, lfy, 'kx')
plt.plot(fx, plfy, 'k-')
plt.annotate('Synth', [0.95*al[1], lfy.mean()])
plt.axis(lal)
plt.ylabel('$\ln (A)$')
plt.title('Amplitude variation with offset')

plt.show()

