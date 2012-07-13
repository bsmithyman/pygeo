#!/usr/bin/env python

import sys
import numpy as np
import struct

import pygeo.segyread as segyread

infile = 'line5-postpick.sgy'
outfile = 'line5-reconstructed.sgy'
scalco = 10.
scalel = 10.

def printnow (text):
  sys.stdout.write(text + '\n')
  sys.stdout.flush()

def overprint (text):
  sys.stdout.write(text + '\r')
  sys.stdout.flush()

printnow('Opening [%s]...'%(infile,))
sf = segyread.SEGYFile(infile, endian='Big')
ns = sf.ns
ntr = sf.ntr
printnow('\tFound %d traces.'%(ntr,))

printnow('Calculating ensemble boundaries...')
sf._calcEnsembles()

# ------------------------------------------------------------------------
# Determine the mapping between SHOTID, gather number and trace number
ngathers = len(sf.ensembles)
ordering = np.argsort(sf.ensembles.values())
shotnums = np.array(sf.ensembles.keys())[ordering]
beginmap = np.array([sf.ensembles[num] for num in shotnums])
printnow('\tFound %d shot gathers.'%(ngathers,))

# ------------------------------------------------------------------------
# Read the receiver geometry and form two arrays containing the results
printnow('Reading receiver geometry information...')
geom = np.array([[trh['gx'],trh['gy']] for trh in sf.trhead], dtype=np.float32)/scalco
x,y = geom.T
cplgeom = x + 1j * y

# ------------------------------------------------------------------------
# Look for unique values along X axis.  This returns three arrays:
#   1. The unique values of X
#   2. Indices in the original array 'x' that contain the values of xu
#   3. Indices in xu for every receiver in x!  This is sort of a big deal.
#cgeom, indexmap, invmap = np.unique(x + 1j * y, return_index=True, return_inverse=True)

seen = {}
cgeom = []
indexmap = []
invmap = []

for (i, item) in enumerate(cplgeom):
  if item in seen:
    invmap.append(seen[item])
  else:
    lcg = len(cgeom)
    seen[item] = lcg
    cgeom.append(item)
    indexmap.append(i)
    invmap.append(lcg)

cgeom = np.array(cgeom)
xu = cgeom.real
yu = cgeom.imag
zu = np.array([sf.trhead[item]['gelev'] for item in indexmap])/scalel

nrecs = len(xu)
printnow('\tFound %d unique receivers.'%(nrecs,))

# ------------------------------------------------------------------------
# Construct a new trace header database, including the zero traces needed
# to accommodate the new receiver numbering
deadtrhead = sf.trhead[0]
deadtrhead['sx'] = 0
deadtrhead['sy'] = 0
deadtrhead['selev'] = 0
deadtrhead['sstat'] = 0
deadtrhead['gx'] = 0
deadtrhead['gy'] = 0
deadtrhead['gelev'] = 0
deadtrhead['gstat'] = 0
deadtrhead['ep'] = 0
deadtrhead['fldr'] = 0
deadtrhead['cdp'] = 0
deadtrhead['cdpt'] = 0
deadtrhead['tracf'] = 0
deadtrhead['tracl'] = 0
deadtrhead['tracr'] = 0
deadtrhead['trid'] = 3

deadtrace = np.zeros((ns,), dtype=np.float32)

bitmaparr = np.zeros((ngathers,nrecs))
livetracenum = 0
totaltracenum = 0

printnow('Remapping headers and writing SEG-Y file...')
with open(outfile, 'w+b') as fp:
  infp = sf._fp
  fp.write(sf.thead.encode('IBM500')[:3200])
  bheadbin = struct.pack('>3L24H', *[sf.bhead[key] for key in segyread.BHEADLIST]) + '\x00' * 340
  fp.write(bheadbin)

  for i in xrange(ngathers):
    overprint('\tSource: %4d\r'%(i+1,))
    start = beginmap[i]
    try:
      stop = beginmap[i+1]
    except IndexError:
      stop = ntr - 1

    templhead = sf.trhead[start]
    localdeadtrhead = deadtrhead.copy()
    localdeadtrhead['fldr'] = templhead['fldr']
    localdeadtrhead['sx'] = templhead['sx']
    localdeadtrhead['sy'] = templhead['sy']
    localdeadtrhead['selev'] = templhead['selev']

    for j in xrange(nrecs):
      totaltracenum += 1
      if (invmap[livetracenum] == j):
        bitmaparr[i][j] += 1

        newtrhead = sf.trhead[livetracenum]
        livetracenum += 1
      else:
        newtrhead = localdeadtrhead.copy()
        newtrhead['gx'] = xu[j]*scalco
        newtrhead['gy'] = yu[j]*scalco
        newtrhead['gelev'] = zu[j]*scalel
        newtrhead['offset'] = np.sqrt(
        			  (xu[j]*scalco-newtrhead['sx'])**2
        			+ (yu[j]*scalco-newtrhead['sy'])**2
				+ (zu[j]*scalel-newtrhead['selev'])**2)

      newtrhead['tracf'] = j + 1
      newtrhead['tracl'] = totaltracenum
      newtrhead['tracr'] = totaltracenum
      newtrhead['ep'] = i + 1

      trheadbin = struct.pack(segyread.TRHEADSTRUCT, *[newtrhead[key] for key in segyread.TRHEADLIST]) + '\x00' * 60
      fp.write(trheadbin)

      if (bitmaparr[i][j] != 0):
        infp.seek(sf._calcDataOffset(livetracenum, ns))
        fp.write(infp.read(ns*4))
      else:
        fp.write(deadtrace)

printnow('Completed writing [%s]!'%(outfile,))
