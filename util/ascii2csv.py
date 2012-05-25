#!/usr/bin/env python

from pygeo.fastio import *
import numpy as np
import scipy as sp
import glob
import multiprocessing
import sys

dims = (301,501,61)
modelfile = 'vel.mod'
modelout = 'velmod.raw'
hcfile = 'num.cell'
hcout = 'numcell.raw'
residout = 'resids.csv'
origfilewc = '*.ascii'
calcfilewc = '*.calcascii'

def writeCSV (filename, items):
  midpoints, otimes, ctimes, errors, offsets = items
  with open(filename, 'w') as fp:
    fp.write('x,y,z,otimes,ctimes,diff,offsets\n')
    for i in xrange(len(otimes)):
      fp.write('%e,%e,%e,%e,%e,%e,%e\n'%(
			midpoints[i][0],
			midpoints[i][1],
			midpoints[i][2],
			otimes[i],
			ctimes[i],
			errors[i],
			offsets[i]))

def syncprint (message):
  sys.stdout.write(message + '\n')
  sys.stdout.flush()


def getTimes (dset):
  return dset['t']


def getMidpoints (dset):
  sc = dset['sc']
  rc = dset['rc']

  return (rc + sc)/2.


def getOffsets (dset):
  sc = dset['sc']
  rc = dset['rc']
  
  return np.sqrt(((rc - sc)**2).sum(axis=1))
  

pool = multiprocessing.Pool()

origfiles = glob.glob(origfilewc)
origfiles.sort()
calcfiles = glob.glob(calcfilewc)
calcfiles.sort()

syncprint('Reading pick files...')
originfo = pool.map(readPicks, origfiles)

syncprint('Reading calc files...')
calcinfo = pool.map(readPicks, calcfiles)

syncprint('\nCalculating:')

syncprint('\ttimes...')
otimes = np.concatenate(map(getTimes, originfo))
ctimes = np.concatenate(map(getTimes, calcinfo))
errors = ctimes - otimes

syncprint('\tmidpoints...')
midpoints = np.concatenate(map(getMidpoints, originfo))

syncprint('\toffsets...')
offsets = np.concatenate(map(getOffsets, originfo))

syncprint('\nWriting output file...')
items = [midpoints, otimes, ctimes, errors, offsets]
writeCSV(residout, items)

syncprint('Done!')
