
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

import mmap
import os
import struct
import numpy as np

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

STRUCT_TRHEAD = '>7L4H8L2H4L46H'
STRUCT_BHEAD = '>3L24H'


class SEGYArray (np.ndarray):
  '''SEG-Y Array'''

  fm = None
  bhead = None
  thead = None
  
  def __new__ (cls, filename):
    filesize = os.path.getsize(filename)

    fp = open(filename, 'r+b')
    fm = mmap.mmap(fp.fileno(), 0)
    fp.close()

    trhead = fm[3600:3840]
    ns = struct.unpack(STRUCT_TRHEAD, trhead[:180])[38]
    
    shape = ((filesize - 3600) / (240 + 4*ns),60+ns)

    self = np.ndarray.__new__(cls, dtype=np.float32, shape=shape, buffer=fm, offset=3600)

    self = self[:,60:]
    self.fm = fm
    self.thead = fm[:3200].replace(' ','\x25').decode('IBM500')
    tbh = struct.unpack(STRUCT_BHEAD, fm[3200:3260])
    self.bhead = {}
    for i, label in enumerate(BHEADLIST):
      self.bhead[label] = tbh[i]

    return self

  def get_trhead(self, trace):
    '''Gets trace headers for a zero-based trace number.'''

    offset = 3600 + trace*self.strides[0]
    tth = struct.unpack(STRUCT_TRHEAD,self.fm[offset:offset+180])
    traceheaders = {}
    for i, label in enumerate(TRHEADLIST):
      traceheaders[label] = tth[i]

    return traceheaders

class SUArray (np.ndarray):
  '''SU Array'''

  fm = None
  bhead = None
  thead = None
  
  def __new__ (cls, filename):
    filesize = os.path.getsize(filename)

    fp = open(filename, 'r+b')
    fm = mmap.mmap(fp.fileno(), 0)
    fp.close()

    trhead = fm[:240]
    ns = struct.unpack(STRUCT_TRHEAD, trhead[:180])[38]
    
    shape = (filesize / (240 + 4*ns),60+ns)

    self = np.ndarray.__new__(cls, dtype=np.float32, shape=shape, buffer=fm)

    self = self[:,60:]
    self.fm = fm

    self.thead = None
    self.bhead = None

    return self

  def get_trhead(self, trace):
    '''Gets trace headers for a zero-based trace number.'''

    offset = 3600 + trace*self.strides[0]
    tth = struct.unpack(STRUCT_TRHEAD,self.fm[offset:offset+180])
    traceheaders = {}
    for i, label in enumerate(TRHEADLIST):
      traceheaders[label] = tth[i]

    return traceheaders

