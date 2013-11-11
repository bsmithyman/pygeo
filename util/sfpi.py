#!/usr/bin/env python
#
# SEG-Y/FAST Pick Importer
# Brendan Smithyman
# bsmithyman@eos.ubc.ca
# March, 2011
#
# sfpi.py
#
# pygeo - a distribution of tools for managing geophysical data
# Copyright (C) 2011, 2012 Brendan Smithyman
#
#                              /---------\
#                              | LICENSE |
#                              \---------/
#
#
# This file is part of pygeo.
#
# pygeo is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# pygeo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with pygeo.  If not, see <http://www.gnu.org/licenses/>.
#
# ----------------------------------------------------------------------
#
#                              /---------\
#                              | READ ME |
#                              \---------/
#
# This program accesses a target SEG-Y file and maps out its headers
# (namely SHOTID and CHANNEL, as specified in the 'headers' dictionary
# below). The SEG-Y file should have contiguous shot gathers, but the
# order of the gathers and the order of traces within the gather does
# not matter. The geometry of the sources and receivers in the SEG-Y
# file is cached (i.e. sx,sy and rx,ry headers for each trace).
#
# The program accesses ASCII-format pick files compatible with FAST
# (First Arrival Seismic Tomography, by Colin Zelt of Rice Unviersity)
# and compares the shot and receiver locations to the cached headers
# from the target SEG-Y file.
#
# Depending on the value of the flag 'insertheaders', this program can
# directly inject the FAST-format pick times into the SEG-Y headers for
# the corresponding trace. Regardless of the value of the flag, a text
# file is written containing the mappings between the SEG-Y SHOTID and
# CHANNEL entries and the FAST picks.
#
# The units of the geometry information are read from the COORD_SCALE
# header of the first trace in the SEG-Y file; however it is necessary
# to choose 'fastunit' to compensate if you are using different units in
# the FAST pick files.
#
# The default distance defined as coincident is 5 m; change 'threshold'
# if you need a different value.

# ----------------------------------------------------------------------
# Module imports
import os
import sys
import numpy as np
import mmap
import struct
import glob
import warnings

# Aliases
ssw = sys.stdout.write
ssf = sys.stdout.flush

# ----------------------------------------------------------------------
# Configuration parameters
pickmask = 'fast/fd*.ascii'
sgfile = sys.argv[1]
outfile = 'convpicks.dat'
insertheaders = False#True

# What is the unit for fast geometry?
fastunit = 1000 # metres or ns

# How close is close?
threshold = 5 # metres

# Mapping between headers and byte-offsets
headers = {'shotid':8,'channel':12,'sx':72,'sy':76,'rx':80,'ry':84,'delay':108}

# ----------------------------------------------------------------------
# Internal constants / variables
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

progchars = ['-','\\','|','/']

# ----------------------------------------------------------------------
# Useful functions

def fastgetshotxy (pickfile):
  f = open(pickfile, 'r')
  line = f.readline()
  f.close()
  return [float(item)*fastunit for item in line.strip().split()[:2]] #x,y

def ssc ():
  sys.stdout.write('                                                                       \r')
  ssf()

def progress ():
  if (progress.count > 3):
    progress.count = 0
  ssw('%s\r'%(progchars[progress.count],))
  ssf()
  progress.count += 1
progress.count = 0

# ----------------------------------------------------------------------
# Main program

print('\nSEG-Y/FAST Pick Importer\nBrendan Smithyman\nbsmithyman@eos.ubc.ca\nMarch, 2011\n')

filesize = os.path.getsize(sgfile)

ssw('  Mapping SEG-Y file... ')
f = open(sgfile, 'r+b')
sgmap = mmap.mmap(f.fileno(),0,mmap.MAP_SHARED)
f.close()
ssw('Done.\r*\n')

if (filesize != sgmap.size()):
  ssw('ERROR: Failed to map entire file.\n')

ssw('  Scanning SEG-Y file and determining parameters... ')
# Decode EBCDIC reel header
textheader = sgmap[:3200].decode('IBM500')

# Decode and unpack binary header
tbh = struct.unpack(STRUCT_BHEAD, sgmap[3200:3260])
bhead = {}
for i, label in enumerate(BHEADLIST):
  bhead[label] = tbh[i]

ns = bhead['hns']

if (ns != struct.unpack('>H',sgmap[3600+114:3600+116])[0]):
  print('ERROR: Number of samples in first trace header conflicts with\n       binary header.')
  exit()

if ((filesize-3600)%(240+ns*4) != 0):
  print('ERROR: File length mismatch.\n       Likely that file is corrupt or has variable-length traces.')
  exit()

ntr = (sgmap.size()-3600)/(240+ns*4)

def calcoffset (trace):
  return 3600 + (240+ns*4)*trace

coordscale = struct.unpack('>h',sgmap[3670:3672])[0]
if (coordscale < 0):
  coordscale = 1./abs(coordscale)

ssw('Done.\r*\n')

# Compile an array of all shotids
shotoff = headers['shotid']
chanoff = headers['channel']
shotids = []
chanids = []
ssw('  Reading SHOTID and CHANNEL headers... \r')
ssf()
for trace in xrange(ntr):
  off = calcoffset(trace)
  shotids.append(struct.unpack('>L',sgmap[off+shotoff:off+shotoff+4])[0])
  chanids.append(struct.unpack('>L',sgmap[off+chanoff:off+chanoff+4])[0])
  progress()

ssw('\r* Reading SHOTID and CHANNEL headers... Done.\n')

shotids = np.array(shotids)
chanids = np.array(chanids)

ssw('  Finding unique SHOTIDs... ')
ssf()
# Find all unique shot ids and an occurence
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    [usids, usidlocs] = np.unique(shotids,return_index=True)
usidorder = np.argsort(usidlocs)
usids = usids[usidorder]
usidlocs = usidlocs[usidorder]
ssw('Done.\r*\n')

sxo = headers['sx']

ssw('  Reading coordinates for each unique SHOTID... ')
ssf()
shotlocations = np.array([[item*coordscale for item in struct.unpack('>2L',sgmap[calcoffset(loc)+sxo:calcoffset(loc)+sxo+8])] for loc in usidlocs])
ssw('Done.\r*\n')
 
rxo = headers['rx']

ssw('  Reading coordinates for each CHANNEL by SHOTID gather... \r')
ssf()
tracedata = []
for sid in usids:
  traces = np.argwhere(shotids == sid)
  tracedata.append(np.array([[loc]+[item*coordscale for item in struct.unpack('>2L',sgmap[calcoffset(loc)+rxo:calcoffset(loc)+rxo+8])] for loc in traces.reshape((len(traces),))]))
  progress()
ssw('* Reading coordinates for each CHANNEL by SHOTID gather... Done.\n')

# Scan all the matching pick files
ssw('  Reading shot coordinates from FAST files...')
ssf()
fastfiles = glob.glob(pickmask)
fastshotinfo = map(fastgetshotxy,fastfiles)
fastorder = np.argsort(fastfiles)
ssw('Done.\r*\n')

outlines = ['SHOTID\tCHANNEL\tTime\n']
fastcorresponds = []
for index in fastorder:
  fsx,fsy = fastshotinfo[index]
  shotidloc = np.argwhere((np.abs(shotlocations[:,0]-fsx) < threshold) * (np.abs(shotlocations[:,1]-fsy) < threshold))
  fastcorresponds.append(shotidloc)
  ssw('\r  Processing %r, which corresponds to SHOTID %d'%(fastfiles[index],usids[shotidloc]))
  ssf()
  f = open(fastfiles[index], 'r')
  lines = f.readlines()
  f.close()

  pickinfo = np.array([[float(item)*fastunit for item in line.strip().split()[:4]] for line in lines[1:]])
  for frx,fry,frz,ptime in pickinfo:
    lokey = np.argwhere((np.abs(tracedata[shotidloc][:,1].reshape((tracedata[shotidloc].shape[0],))-frx) < threshold) * (np.abs(tracedata[shotidloc][:,2].reshape((tracedata[shotidloc].shape[0],))-fry) < threshold))
    if (lokey):
      trace = int(tracedata[shotidloc][lokey,0][0,0])
      outlines.append('%d\t%d\t%f\n'%(shotids[trace],chanids[trace],ptime))
      if (insertheaders):
        delayloc = calcoffset(trace) + headers['delay']
        sgmap[delayloc:delayloc+2] = struct.pack('>H',int(ptime))

ssw('\n  Processed picks...\n')
ssf()
if (insertheaders):
  ssw('(flushing buffers)...\n')
  sgmap.flush()
  sgmap.close()
  ssw('                 ...and inserted them into the trace headers.\n')

ssw('  Opening %s...\n'%(outfile,))
ssf()
f = open(outfile,'w')
f.writelines(outlines)
f.close()
ssw('  Wrote %s.\n'%(outfile,))
ssw('\nComplete!\n')
