#!/usr/bin/ipython2 -i

# pygeo - a distribution of tools for managing geophysical data
# Copyright (C) 2011, 2012, 2013 Brendan Smithyman

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

import os
import sys
import warnings

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from pygeo.segyread import SEGYFile
import pygeo.analysis as analysis

warnings.filterwarnings('ignore')

infile = 'bp_data2_data.sgy'
tempfile = 'tmp/traces.dat'
watervelocity = 1482.
keepdata = False

dbase = {}

# ------------------------------------------------------------------------
# Functions

def writenow (text):
  sys.stdout.write(text + '\n')
  sys.stdout.flush()

def loadfunc ():

  sf = SEGYFile(infile)
  sf._calcEnsembles()
  ntr = sf.ntr
  nsrc = len(sf.ensembles)
  nrec = ntr / nsrc
  nsamp = sf.ns
  dt = sf.bhead['hdt'] / 1.e6

  trh0 = sf.trhead[0]

  scalco = trh0['scalco']
  if (scalco < 0):
    scalco = -1. / scalco

  scalel = trh0['scalel']
  if (scalel < 0):
    scalel = -1. / scalel

  datadims = (nsrc,nrec,nsamp)
  headerdims = (nsrc,nrec)
  flatdatadims = (ntr,nsamp)

  params = {
    'ntr': ntr,
    'nsrc': nsrc,
    'nrec': nrec,
    'nsamp': nsamp,
    'dt': dt,
    'scalco': scalco,
    'scalel': scalel,
    'datadims': datadims,
    'headerdims': headerdims,
    'flatdatadims': flatdatadims,
  } 

  dbase['params'] = params

  sx = np.array([float(abs(trh['sx'])) for trh in sf.trhead]).reshape(headerdims)
  gx = np.array([float(abs(trh['gx'])) for trh in sf.trhead]).reshape(headerdims)
  offset = abs(sx - gx)
  cdp = (sx + gx)/2.

  water = offset / watervelocity
  wmute = water / dt

  headers = {
    'sx': sx,
    'gx': gx,
    'offset': offset,
    'cdp': cdp,
    'water': water,
    'wmute': wmute,
  }

  dbase['headers'] = headers

  if (os.path.exists(tempfile) and os.path.getmtime(tempfile) > os.path.getmtime(infile)):
    tracebuf = np.memmap(tempfile, dtype=np.float32, shape=datadims, mode='r+')
  else:
    tracebuf = np.memmap(tempfile, dtype=np.float32, shape=datadims, mode='w+')
    tracebuf[:] = sf[:].reshape(datadims)

  dbase['data_orig'] = tracebuf
  dbase['data_current'] = tracebuf

def process_reset ():
  writenow('PROCESS: Resetting all data arrays.')
  nukekeys = []
  for key in dbase:
    if (key.find('data_') == 0 and key != 'data_orig'):
      nukekeys.append(key)

  for key in nukekeys:
    del dbase[key]

  dbase['data_current'] = dbase['data_orig']

def process_mutebackscatter (shift = 0., inplace = False):
  writenow('PROCESS: Muting reflection section after direct wave.')
  writenow('         \tshift = %f s'%shift)
  params = dbase['params']
  headers = dbase['headers']
  data = dbase['data_current']
  dt = params['dt']
  nsamp = params['nsamp']
  tmax = dt*nsamp
  nsrc = params['nsrc']
  nrec = params['nrec']
  wmute = headers['wmute']

  if (not inplace):
    data = data.copy()

  for i in xrange(nsrc):
    for j in xrange(nrec):

      data[i,j,max(0,int(wmute[i,j]+shift/dt)):] = 0.

  if (keepdata):
    dbase['data_mutebackscatter'] = data
  dbase['data_current'] = data

def process_energyratio (windowsize = 0.1):
  writenow('PROCESS: Computing STA/LTA for traces.')
  writenow('         \twindowsize = %f s'%windowsize)
  params = dbase['params']
  data = dbase['data_current'].reshape(params['flatdatadims'])
  dt = params['dt']

  windowsamps = int(windowsize/dt)
  erdata = analysis.energyratio(data, windowsamps).reshape(params['datadims'])

  if (keepdata):
    dbase['data_energyratio'] = erdata
  dbase['data_current'] = erdata

def process_normalize ():
  writenow('PROCESS: Normalizing trace amplitudes.')

  params = dbase['params']
  data = dbase['data_current'].reshape(params['flatdatadims'])
  
  normdata = analysis.tracenormalize(data).reshape(params['datadims'])

  if (keepdata):
    dbase['data_normalize'] = normdata
  dbase['data_current'] = normdata

def process_agc (windowsize = 0.1):
  writenow('PROCESS: Applying Automatic Gain Control.')
  writenow('         \twindowsize = %f s'%windowsize)
  params = dbase['params']
  data = dbase['data_current'].reshape(params['flatdatadims'])
  dt = params['dt']

  windowsamps = int(windowsize/dt)
  agcdata = analysis.agc(data, windowsamps).reshape(params['datadims'])

  if (keepdata):
    dbase['data_energyratio'] = agcdata
  dbase['data_current'] = agcdata

def display_gather (shotid):
  writenow('DISPLAY: Showing shot gather with header(s) overlaid.')
  writenow('         \tSHOTID = %d'%shotid)
  params = dbase['params']
  headers = dbase['headers']
  data = dbase['data_current']

  tmax = params['dt']*params['nsamp']
  nrec = params['nrec']
  water = headers['water']

  fig = plt.figure()
  ax = fig.add_subplot(1,1,1, aspect='auto')

  gather = data[shotid,:]
  bound = max(abs(gather.min()),abs(gather.max()))

  ax.imshow(gather.T, extent=[0, nrec, tmax, 0], aspect='auto', cmap=cm.bwr, vmin=-bound, vmax=bound)
  ax.plot(np.arange(0.5, nrec+0.5), water[shotid,:], 'k-', label='Direct Wave')

  for key in headers:
    if (key.find('picks_') == 0):
      pickname = key[6:]
      ax.plot(np.arange(0.5, nrec+0.5), headers[key][shotid,:], label=pickname)

  plt.axis([0.5, nrec+0.5, tmax, 0])
  plt.legend()

  plt.show()
  
def pick_argwhere (shift = 0.):
  writenow('PICKING: Using location of first non-zero entry in each trace.')
  writenow('         \tshift = %f s'%shift)
  params = dbase['params']
  headers = dbase['headers']
  data = dbase['data_current'].reshape(params['flatdatadims'])

  dt = params['dt']

  picks = []
  for trace in data:
    try:
      loc = float(np.argwhere(trace != 0)[0])
    except:
      loc = 0.

    picks.append(dt*loc + shift)

  picks = np.array(picks).reshape(params['headerdims'])

  headers['picks_argwhere'] = picks

def headers_masknearoffset (inkey, offset):
  writenow('HEADERS: Masking header at near offsets.')
  writenow('         \theader = %s, offset = %f m'%(inkey,offset))
  headers = dbase['headers']

  mask = headers['offset'] < offset
  headers[inkey + '_masked'] = np.ma.MaskedArray(headers[inkey], mask)

# ------------------------------------------------------------------------
# Main Program

loadfunc()

params = dbase['params']
headers = dbase['headers']

#pick_argwhere(0.5)
#headers_masknearoffset('picks_argwhere', 5000)
#plt.imshow(headers['picks_argwhere_masked'], aspect='auto')
#plt.show()
