#!/usr/bin/env python

import numpy as np
from pygeo.segyread import SEGYFile
from pygeo.analysis import energyratio
from pygeo.fullpy import readini
import glob
import sys

DEFAULT_WLEN = 10
TUNIT = 1e-3 # milliseconds
# TUNIT = 1e-6 # seconds

try:
    infile = sys.argv[1]
except IndexError:
    raise Exception('No input file provided!')

wlen = DEFAULT_WLEN

if len(sys.argv) >= 3:
    wlen = int(sys.argv[2])

shift = -wlen

if len(sys.argv) >= 4:
    shift = -int(sys.argv[3])

inifile = glob.glob('*.ini')[0]
projnm = inifile.split('.')[0]

ini = readini(inifile)

ns = ini['ns']
nr = ini['nr']

sf = SEGYFile(infile)
dt = sf.bhead['hdt'] * TUNIT
picks = (np.argmax(energyratio(sf[:], wlen), axis=1) + shift)
picks.shape = (ns, nr)

picks = picks * dt

with open('%s.picks'%(projnm,), 'w') as fp:
    fp.writelines(['%d %d %f\n'%(isrc+1, irec+1, picks[isrc, irec]) for isrc in xrange(ns) for irec in xrange(nr)])


