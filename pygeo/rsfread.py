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

import os
import numpy as np
import re
import sys
import mmap

PARSE_EXPR = '[\n\t ]*(?P<param>[a-zA-Z0-9]+)=[\"\\\']?(?P<value>[a-zA-Z0-9\.\@]+)[\"\\\']?[\n\t ]*'
HEADER_END = '\x0c\x0c\x04'

parse_finder = re.compile(PARSE_EXPR).findall
getopts = lambda instr: dict(parse_finder(instr))

dtype_mapping = {
	'native': np.float32,
	'native_float': np.float32,
	'float': np.float32,
	'complex': np.complex64,
	'double': np.float64,
}

maxhlen = 1024


class RSFFile (np.memmap):
  '''
  Provides read access to a RSF dataset (headers and data).
  '''

  headfile = None
  datafile = None
  verbose = False
  usemmap = True

  header = None

  # --------------------------------------------------------------------

  def _maybePrint (self, text):
    if self.verbose:
      sys.stdout.write('%s\n'%(text,))
      sys.stdout.flush()

  # --------------------------------------------------------------------

  def __new__ (cls, filename, verbose = None, usemmap = None, mode = 'r+'):

    if (verbose is not None):
      self.verbose = verbose

    if (usemmap is not None):
      self.usemmap = usemmap

    self._maybePrint('Opening RSF header file...')

    with open(filename, 'r') as fp:
      headbin = fp.read(maxhlen)
      headend = headbin.find(HEADER_END)

      if (headend != -1):
        offset = headend + len(HEADER_END)
        header = getopts(headbin[:headend])
      else:
        offset = 0
        header = getopts(headbin)

    if (offset == 0):
      try:
        datafile = header['in']
        self._maybePrint('... got data file "%s".'%(self.datafile,))
      except:
        print('RSF header does not specify data file!')
        raise

      if (not os.path.isfile(datafile)):
        print('RSF data file "%s" does not exist!'%(datafile,))

    else:
      datafile = filename
      self._maybePrint('... data in same file at offset %d.'%(offset,))

    datalen = os.path.getsize(datafile) - offset
    fp = open(datafile, 'r+b')

    if (usemmap):
      try:
        self._maybePrint('Trying to create memory map...')
        _fp = mmap.mmap(fp.fileno(), 0)
        self._maybePrint('Success. Using memory-mapped I/O.\n')
      except:
        _fp = fp
        usemmap = False
        #self._maybePrint('Memory map failed; using conventional I/O.\n')
    else:
      _fp = fp

    dimlist = []
    for key in header:
    
      if (len(key) == 2 and key[0] == 'n'):
        header[key] = int(header[key])
        dimlist.append(key)

      elif (len(key) == 2 and (key[0] == 'd' or key[0] == 'o')):
        header[key] = float(header[key])

      elif (key == 'esize'):
        header[key] == int(header[key])
        bytesize = header[key]

    dimlist.sort(reverse=True)
    dims = tuple([header[dim] for dim in dimlist])

    self = np.memmap.__new__(cls, datafile, dtype=dtype_mapping[header['format']], offset=offset, mode=mode)
    self.shape = dims
    self.header = header
    self.datafile = datafile
    self.headfile = filename

    return self

