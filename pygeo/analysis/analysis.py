
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

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

interpscalco = lambda scalco: float(-1./scalco) if (scalco < 0) else float(scalco)

calc2Doffset = lambda trh: np.sqrt((trh['gx'] - trh['sx'])**2 + (trh['gy'] - trh['sy'])**2) * interpscalco(trh['scalco'])

calc3Doffset = lambda trh: np.sqrt(((trh['gx'] - trh['sx'])**2 + (trh['gy'] - trh['sy'])**2) * interpscalco(trh['scalco'])**2 + ((trh['gelev'] - trh['selev']) * interpscalco(trh['scalel']))**2)

getlive = lambda offsets: [i for i in xrange(len(offsets)) if offsets[i] != 0]

timereduce = lambda offsets, redvel, shift: [float(offset) / redvel + shift for offset in offsets]

def getpicks (sf, header='delrt'):
  picks = np.array([(float(trh[header]) if trh[header] < 60000 else 0.) for trh in sf.trhead])
  return picks

def calcoffset (sf, in3D=False):
  norm = (calc3Doffset if in3D else calc2Doffset)
  return np.array([norm(trh) for trh in sf.trhead])

def clipsign (value, clip):
  clipthese = abs(value) > clip
  return value * ~clipthese + np.sign(value)*clip*clipthese

def wiggle (traces,skipt=1,scale=1.,lwidth=.1,offsets=None,redvel=0.,tshift=0.,sampr=1.,clip=10.,color='black',fill=True,line=True):

  if (offsets is not None):
    shifts = timereduce(offsets, redvel, tshift)
  else:
    shifts = np.zeros((ntr,))

  ns = traces.shape[1]
  ntr = traces.shape[0]
  t = np.arange(ns)*sampr

  for i in range(0, ntr, skipt):
    trace = traces[i].copy()
    trace[0] = 0
    trace[-1] = 0

    if (line):
      plt.plot(i + clipsign(trace / scale, clip), t - shifts[i], color=color, linewidth=lwidth)
    if (fill):
      for j in range(ns):
        if (trace[j] < 0):
          trace[j] = 0
      plt.fill(i + clipsign(trace / scale, clip), t - shifts[i], color=color, linewidth=0)
  plt.grid(True)

def tracenormalize (traces):
  dims = traces.shape
  nsamp = dims[-1]
  if (dims[:-1] > 1):
    ntr = reduce(lambda x,y: x*y, dims[:-1])
  else:
    ntr = dims[0]

  result = np.array([trace/max(abs(trace.max()),abs(trace.min())) for trace in traces.reshape(ntr,nsamp)]).reshape(dims)
  return result

def agc (traces, windowlen):
  dims = traces.shape
  nsamp = dims[-1]
  if (dims[:-1] > 1):
    ntr = reduce(lambda x,y: x*y, dims[:-1])
  else:
    ntr = dims[0]

  tracetemp = traces.reshape(ntr,nsamp)
  result = np.zeros_like(tracetemp)

  windowshift = windowlen/2

  for i in xrange(ntr):
    windowsum = (tracetemp[i,:windowshift]**2).sum()
    for j in xrange(windowshift,nsamp):

      result[i,j] = windowlen * tracetemp[i,j] / windowsum

      try:
        addenergy = tracetemp[i,j+windowshift]**2
      except IndexError:
        addenergy = 0.

      try:
        subenergy = tracetemp[i,j-windowshift]**2
      except:
        subenergy = 0.

      windowsum += addenergy
      windowsum -= subenergy

  return result

