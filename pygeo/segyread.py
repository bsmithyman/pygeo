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

import os.path
import mmap
import struct
import numpy as np
import scipy as sp
import sys

MEGABYTE = 1048576


BHEADLIST = ['jobid','lino','reno','ntrpr','nart','hdt','dto','hns','nso',
             'format','fold','tsort','vscode','hsfs','hsfe','hslen','hstyp',
             'schn','hstas','hstae','htatyp','hcorr','bgrcv','rcvm','mfeet',
             'polyt','vpol']

TRHEADLIST = ['tracl','tracr','fldr','tracf','ep','cdp','cdpt','trid','nvs',
             'nhs','duse','offset','gelev','selev','sdepth','gdel','sdel',
             'swdep','gwdep','scalel','scalco','sx','sy','gx','gy','counit',
             'wevel','swevel','sut','gut','sstat','gstat','tstat','laga',
             'lagb','delrt','muts','mute','ns','dt','gain','igc','igi','corr',
             'sfs','sfe','slen','styp','stas','stae','tatyp','afilf','afils',
             'nofilf','nofils','lcf','hcf','lcs','hcs','year','day','hour','minute','sec',
             'timbas','trwf','grnors','grnofr','grnlof','gaps','otrav']

MAJORHEADERS = [1,2,3,7,38,39]
 

class SEGYFile:
  '''
  Provides read access to a SEG-Y dataset (headers and data).
  '''

  mmaplimit = 1000*MEGABYTE # Limit memory mapped i/o to files under 1000M

  filename = None
  verbose = True
  majorheadersonly = True
  isSU = False
  endian = 'Auto'

  samplen = 4

  mendian = None
  usemmap = False
  thead = None
  bhead = None
  trhead = None
  ntr = None
  ensembles = None
  initialized = False

  # --------------------------------------------------------------------

  # Written by Robert Kern on the SciPy-user mailing list.
  def _ibm2ieee (self, ibm):
    """ Converts an IBM floating point number into IEEE format. """

    sign = ibm >> 31 & 0x01

    exponent = ibm >> 24 & 0x7f

    mantissa = ibm & 0x00ffffff
    mantissa = (mantissa * 1.0) / pow(2, 24)

    ieee = (1 - 2 * sign) * mantissa * pow(16.0, exponent - 64)

    return ieee

  # --------------------------------------------------------------------

  def _detect_machine_endian (self):
    '''
    Detects native (machine) endian.
    '''

    if struct.pack('h', 1) == '\x01\x00':
      endian = 'Little'
    else:
      endian = 'Big'

    return endian

  # --------------------------------------------------------------------

  def _isInitialized (self):
    return self.initialized

  # --------------------------------------------------------------------

  def _maybePrint (self, text):
    if self.verbose:
      sys.stdout.write('%s\n'%(text,))

  # --------------------------------------------------------------------

  def fOpen (self, filename):
    if self.usemmap:
      return self._map
    else:
      return open(filename, 'rb')

  def fClose (self, fp):
    if not self.usemmap:
      fp.close()
    return

  # --------------------------------------------------------------------

  def _readHeaders (self):

    '''
    Read SEG-Y headers.

    Returns 3 elements:

      1. Text header:   Returned as ASCII text, converted from IBM500 EBCIDIC.
      2. Block header:  Returned as a dictionary of values, named per SU style.
      3. Trace headers: Returned as a list of dictionaries, in trace order.
                        Dictionaries of values are named per SU style.
    '''

    self._maybePrint('Reading SEG-Y headers...')
    fp = self.fOpen(self.filename)

    if (not self.isSU):
      textheader = fp.read(3200).replace(' ','\x25').decode('IBM500')
      blockheader = fp.read(400)

      blockheader = struct.unpack('>3L24H',blockheader[:60])
      bhead = {}

      for i, label in enumerate(BHEADLIST):
        bhead[label] = blockheader[i]

    else:
      textheader = None
      bhead = None

    traceheaders = []

    try:
      while True:
        traceheader = fp.read(240)
        traceheader = struct.unpack('>7L4H8L2H4L46H',traceheader[:180])
        tracehead = {}

        if (self.majorheadersonly):
          for i in MAJORHEADERS:
            tracehead[TRHEADLIST[i]] = traceheader[i]
        else:
          for i, label in enumerate(TRHEADLIST):
            tracehead[label] = traceheader[i]

        traceheaders.append(tracehead)
        fp.seek(4*tracehead['ns'],1)

    except:
      self.fClose(fp)

    self._maybePrint('Read SEG-Y headers.\n\t%d traces present.\n' % (len(traceheaders)))

    [self.thead, self.bhead, self.trhead, self.ntr] = [textheader, bhead, traceheaders, len(traceheaders)]
    return

  # --------------------------------------------------------------------

  def _calcOffset (self, trace, ns):
    '''
    Calculates the byte offset of the beginning of the data portion of a
    seismic trace, given the trace number and the number of samples per.
    '''

    if (not self.isSU):
      return 3200 + 400 + 240 + (ns*self.samplen + 240)*(trace-1)
    else:
      return 240 + (ns*4 + 240)*(trace-1)

  # --------------------------------------------------------------------

  def _detectEndian (self):
    if (self.endian != 'Auto'):
      self._maybePrint('%s endian specified... Not autodetecting.'%(self.endian,))
      if (self.endian != self.mendian):
        self._maybePrint('%s endian != %s endian, therefore Foreign.'%(self.endian,self.mendian))
        self.endian = 'Foreign'
    else:
      self._maybePrint('Auto endian specified... Trying to autodetect data endianness.')
      for i in xrange(1, self.ntr+1):
        locar = self.readTraces(i)
        if ((not abs(locar).sum() == 0.) and (not np.isnan(locar.mean()))):
          nexp = abs(np.frexp(locar.mean())[1])
          locar = locar.newbyteorder()
          fexp = abs(np.frexp(locar.mean())[1])
          if (fexp > nexp):
            self.endian = 'Native'
          else:
            self.endian = 'Foreign'
          self._maybePrint('Scanned %d trace(s). Endian appears to be %s.'%(i, self.endian))
          break

      if (self.endian == 'Foreign'):
        self._maybePrint('Will attempt to convert to %s endian when traces are read.\n'%(self.mendian,))
      elif (self.endian == 'Auto'):
        self._maybePrint('Couldn\'t find any non-zero traces to test!\nAssuming Native endian.\n')


  # --------------------------------------------------------------------

  def readTraces (self, traces=None):
    '''
    Returns trace data as a list of numpy arrays (i.e. non-adjacent trace
    numbers are allowed). Requires that traces be fixed length.
    '''

    ns = self.trhead[0]['ns']
    fp = self.fOpen(self.filename)

    result = []

    if (traces == None):
      traces = range(1, self.ntr+1)

    if not np.iterable(traces):
      traces = [traces]

    # Handles SU format and IEEE floating point
    if (self.isSU or self.bhead['format'] == 5):
      for trace in traces:
        fp.seek(self._calcOffset(trace, ns))
        tracetemp = fp.read(ns*4)
        result.append(np.array(struct.unpack('>%df'%(ns,), tracetemp), dtype=np.float32))

    # Handles everything else
    else:
      if (self._isInitialized()):
        self._maybePrint('FORMAT == %d'%(self.bhead['format'],))

      # format == 1: IBM Floating Point
      if (self.bhead['format'] == 1):
        if (self._isInitialized()):
          self._maybePrint('             ...converting from IBM floating point.\n')
        for trace in traces:
          fp.seek(self._calcOffset(trace, ns))
          tracetemp = struct.pack('%df'%(ns,),*[self._ibm2ieee(item) for item in struct.unpack('>%dL'%(ns,),fp.read(ns*4))])
          result.append(np.array(struct.unpack('>%df'%(ns,), tracetemp), dtype=np.float32))

      elif (self.bhead['format'] == 2):
        if (self._isInitialized()):
          self._maybePrint('             ...reading from 32-bit fixed point.\n')
        for trace in traces:
          fp.seek(self._calcOffset(trace, ns))
          result.append(np.array(struct.unpack('>%dl'%(ns,),fp.read(ns*4)), dtype=np.int32))

      elif (self.bhead['format'] == 3):
        if (self._isInitialized()):
          self._maybePrint('             ...reading from 16-bit fixed point.\n')
        for trace in traces:
          fp.seek(self._calcOffset(trace, ns))
          result.append(np.array(struct.unpack('>%dh'%(ns,),fp.read(ns*2)), dtype=np.int32))

      elif (self.bhead['format'] == 8):
        if (self._isInitialized()):
          self._maybePrint('             ...reading from 8-bit fixed point.\n')
        for trace in traces:
          fp.seek(self._calcOffset(trace, ns))
          result.append(np.array(struct.unpack('>%db'%(ns,),fp.read(ns)), dtype=np.int32))

      elif (self.bhead['format'] == 4):
        if (self._isInitialized()):
          self._maybePrint('             ...converting from 32-bit fixed point w/ gain.\n')
        for trace in traces:
          fp.seek(self._calcOffset(trace, ns))
          tracemantissa = np.array(struct.unpack('>%s'%(ns*'xxh',), fp.read(ns)), dtype=np.float32)
          traceexponent = np.array(struct.unpack('>%s'%(ns*'xbxx',), fp.read(ns)), dtype=np.byte)
          result.append(tracemantissa**traceexponent)

    self.fClose(fp)

    if (self.endian == 'Foreign'):
      return np.array(result, dtype=np.float32).byteswap()
    else:
      return np.array(result, dtype=np.float32)

  # --------------------------------------------------------------------

  def lockTraces (self, traces):
    '''
    Returns trace data as a list of memory-mapped array-like objects.
    Otherwise similar to readTraces.
    '''

    ns = self.trhead[0]['ns']

    result = []

    if not np.iterable(traces):
      traces = [traces]

    for trace in traces:
      tracetemp = np.memmap(self.filename, 'float32', 'c', 
                         offset=self._calcOffset(trace, ns)/4, shape=(ns,))
      result.append(tracetemp)

    return result

  # --------------------------------------------------------------------

  def findTraces (self, key, kmin, kmax):
    if not self.trhead[0].has_key(key):
      raise Invalid_Trace_Header

    validtraces = []

    for i,trace in enumerate(self.trhead):
      if (trace[key] <= kmax and trace[key] >= kmin):
        validtraces.append(i+1)

    return validtraces

  # --------------------------------------------------------------------

  def _calcEnsembles (self):
    self.ensembles = {}

    self._maybePrint('Scanning ensembles...')
    for i, headers in enumerate(self.trhead):

      try:
        self.ensembles.keys().index(headers['fldr'])

      except ValueError:
        self.ensembles[headers['fldr']] = i

    self._maybePrint('Complete. Found %d ensemble(s).\n'%(len(self.ensembles),))
       
  # --------------------------------------------------------------------

  def __init__ (self, filename, verbose = verbose, majorheadersonly = majorheadersonly, isSU = isSU, endian = endian):
    self.filename = filename
    self.verbose = verbose
    self.majorheadersonly = majorheadersonly
    self.isSU = isSU
    self.endian = endian

    self._maybePrint('Detecting machine endianness...')
    self.mendian = self._detect_machine_endian()
    self._maybePrint('%s.\n'%(self.mendian,))

    if os.path.getsize(filename) < self.mmaplimit:
      self.usemmap = True
      self._fp = open(self.filename, 'r+b')
      try:
        self._maybePrint('Trying to create memory map...')
        self._map = mmap.mmap(self._fp.fileno(), 0)
        self._maybePrint('Success. Using memory-mapped I/O.\n')
      except:
        self._fp.close()
        del self._fp
        self.usemmap = False
        self._maybePrint('Memory map failed; using conventional I/O.\n')
    else:
      self.usemmap = False
      self._maybePrint('File exceeds %d MB; using conventional I/O.\n'%(self.mmaplimit / MEGABYTE,))

    # Get header information from file
    self._readHeaders()

    # Determine length of each sample from FORMAT code
    self._getSamplen()

    # Attempt to find shot-record boundaries
    self._calcEnsembles()

    # Autodetect data endian
    self._detectEndian()

    # Confirm that the SEGYFile object has been initialized
    self.initialized = True
  
  # --------------------------------------------------------------------

  def __del__ (self):
    if self.usemmap:
      self._map.close()
      self._fp.close()

  # --------------------------------------------------------------------

  def _getSamplen (self):
    if (self.isSU):
      self.samplen = 4
      return

    if (self.bhead['format'] == 3):
      self.samplen = 2
    elif (self.bhead['format'] == 8):
      self.samplen = 1
    else:
      self.samplen = 4

  # --------------------------------------------------------------------

  def sNormalize (self, traces):
    if not np.iterable(traces):
      traces = [traces]

    self._maybePrint('Normalizing each trace to unit amplitude.\n')

    return np.array([2*(trace - trace.min())/max(trace.max(),abs(trace.min())) - 1 for trace in traces])

  # --------------------------------------------------------------------

  def writeFlat (self, outfilename):

    ntraces = len(self.trhead)

    ns = self.trhead[0]['ns']
    fp = self.fOpen(self.filename)

    fp_out = open(outfilename, "w")

    for trace in range(1, ntraces+1):
      fp.seek(self._calcOffset(trace,ns))
      fp_out.write(fp.read(ns*4))

    fp_out.close()
    self.fClose(fp)

  # --------------------------------------------------------------------

  def writeSEGY (self, outfilename, traces, headers=None):

    if (headers == None):
      thead=self.thead
      bhead=self.bhead
      trhead=self.trhead
    else:
      [thead, bhead, trhead] = headers

    ntraces = len(traces)

    ns = trhead[0]['ns']

    fp = open(outfilename, 'w+b')

    fp.write(thead.encode('IBM500')[:3200])
    
    bheadbin = struct.pack('>3L24H', *[bhead[key] for key in BHEADLIST]) + '\x00' * 340

    fp.write(bheadbin)

    for i, trace in enumerate(traces):
      trheadbin = struct.pack('>7L4H8L2H4L46H', *[trhead[i][key] for key in TRHEADLIST]) + '\x00' * 60
      fp.write(trheadbin)
      tracetemp = struct.pack('>%df'%(ns,), *list(trace))
      fp.write(tracetemp)

    fp.close()

  # --------------------------------------------------------------------

  def writeSU (self, outfilename, traces, trhead=None):

    if (trhead == None):
      trhead=self.trhead

    ntraces = len(traces)

    ns = trhead[0]['ns']

    fp = open(outfilename, 'w+b')

    for i, trace in enumerate(traces):
      trheadbin = struct.pack('>7L4H8L2H4L46H', *[trhead[i][key] for key in TRHEADLIST]) + '\x00' * 60
      fp.write(trheadbin)
      tracetemp = struct.pack('>%df'%(ns,), *list(trace))
      fp.write(tracetemp)

    fp.close()
