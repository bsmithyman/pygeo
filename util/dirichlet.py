#!/usr/bin/env python

import os
import glob
import numpy as np
import matplotlib.pylab as pl
from pygeo.fullpy import readini

# ------------------------------------------------------------------------
# Settings

AUTHORSHIP = 'Brendan Smithyman, June 2012'
VERSION = '%prog v1.0\n'
DESCRIPTION = '''
Computes Dirichlet kernels corresponding to the choices of y-wavenumber used
in an OMEGA configuration file (*.ini).
'''

def_nky = 40
def_vmin = 1500.
def_omega = 8.
def_samps = 1000
def_ymax = 100000

# ------------------------------------------------------------------------
# Functions

def gauleg (m):
  n = (m-1)*2 + 1
  x = np.zeros((m,))
  w = np.zeros((m,))
  eps = 1E-5

  for i in xrange(1,m+1):
    z = np.cos(np.pi*(i-0.25)/(n+0.5))

    while (True):
      p1 = 1.0
      p2 = 0.0
      for j in xrange(1,n+1):
        p3 = p2
        p2 = p1
        p1 = ((2.*j - 1.)*z*p2 - (j-1)*p3)/j

      pp = n*(z*p1 - p2)/(z**2 - 1.)
      z1 = z
      z = z1 - p1/pp
      if (abs(z-z1) <= eps):
        break

    x[m-i] = z
    w[m-i] = 2./((1. - z*z)*pp**2)

  return x, w

def kernelfunc (k, t):
  return np.cos(k * t)

def dirichlet (y, kys, weights=None):
  if (weights is None):
    weights = np.ones_like(kys) / len(kys)

  outseries = 0.5*kernelfunc(kys[0], y)*weights[0]
  for i, ky in enumerate(kys[1:]):
    outseries += kernelfunc(ky, y) * weights[i+1]

  return outseries

def comparison (y, nky, vmin, omega):
  kc = omega/vmin

  lin = np.linspace(0, kc, nky)
  linkernel = dirichlet(y, lin)

  gaux, gauw = gauleg(nky)
  gaux *= kc
  gaukernel = dirichlet(y, gaux, gauw)

  return linkernel, gaukernel

# ------------------------------------------------------------------------
# Parameter parsing

from optparse import OptionParser

usage = '%prog [options] [projectname]'

parser = OptionParser(usage       = usage,
                      version     = VERSION,
                      description = DESCRIPTION)

parser.add_option('-v', '--verbose', action='store_true', dest='verbose',
		help='display additional information')

parser.add_option('-m', '--vmin', action='store', dest='vmin',
		help='minimum velocity for k_c calculation')

parser.add_option('-n', '--nky', action='store', dest='nky',
		help='number of wavenumber components')

parser.set_defaults(
			verbose		= False,
			vmin		= None,
			nky		= None,
)

(options, args) = parser.parse_args()

if (len(args) == 0):
  globresult = glob.glob('*.ini')
  if (len(globresult) == 1):
    projnm = os.path.splitext(globresult[0])[0]
  else:
    parser.error('Please specify a project name!')
else:
  projnm = args[0]

try:
  inidict = readini(projnm + '.ini')
except:
  parser.error('Project ini file %s.ini does not exist!'%(projnm,))

# Decide what the maximum y-difference is and scale accordingly
itemlist = ['srcs','recs','geos']
for item in itemlist:
  if (inidict[item].shape[0] == 0):
    itemlist.remove(item)
  
ymax = max((inidict[item][:,1].max() for item in itemlist))
ymin = min((inidict[item][:,1].min() for item in itemlist))
ydiff = ymax - ymin

# Get existing ky information
kys = inidict['kys']
method = inidict['method']
nky = inidict['nky']
vmin = inidict['vmin']
freqs = inidict['freqs']

if (options.vmin is not None):
  vmin = float(options.vmin)

if (options.nky is not None):
  nky = float(options.nky)

# Define basis for plot
y = np.linspace(0, def_ymax, def_samps)

minomega = freqs.min()
maxomega = freqs.max()

pl.figure()
ex = dirichlet(y, kys)
lk1, gk1 = comparison(y, nky, vmin, minomega)
lk2, gk2 = comparison(y, nky, vmin, maxomega)

pl.subplot(2,1,1)
pl.plot(y, ex, 'r-', label='Current')
pl.plot(y, lk1, 'b-', label='Linear')
pl.plot(y, gk1, 'g-', label='Gauss-Legendre')
pl.ylabel('Amplitude')
pl.xlabel('Cross-line Distance (m)')
pl.title('Minimum $\omega$: %f'%(minomega,))
pl.legend()

pl.subplot(2,1,2)
pl.plot(y, ex, 'r-', label='Current')
pl.plot(y, lk2, 'b-', label='Linear')
pl.plot(y, gk2, 'g-', label='Gauss-Legendre')
pl.ylabel('Amplitude')
pl.xlabel('Cross-line Distance (m)')
pl.title('Maximum $\omega$: %f'%(maxomega,))
pl.legend()

pl.show()
