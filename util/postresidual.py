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
import sys
import glob
import warnings
import numpy as np
from optparse import OptionParser
from pygeo.segyread import SEGYFile

# ------------------------------------------------------------------------
# Settings

AUTHORSHIP = 'Brendan Smithyman'
VERSION = '%prog v1.0\n'
DESCRIPTION = 'Calculates data residuals after inversion, consistent with FULLWV output, but based on final utobs and utest files. Presently does not implement source/receiver weighting.'
USAGE = '%prog [projnm]'
EPILOG = '''
Norms:
   0: conventional amp + phase
   1:              amplitude normalized phase only
   2:              amplitude preserved phase only
  -1:  logarithmic phase only
  -2:              phase + amplitude
  -3:              amplitude only
'''

warnings.simplefilter('ignore', RuntimeWarning)

mintime = 0.01
maxtime = 1.0E10

gen_obsfilekey = '*.utobs*'
gen_estfilekey = '*.utest*'
proto_obsfilekey = '%s.udiff*'
proto_estfilekey = '%s.utest*'

validnorms = [-3,-2,-1,0,1,2]
picknorms = [0,2]
nantag = float('339999995214436424907732413799364296704.0000')
refloc = '>TRCGEOMETRY\n'
refoffset = 2
tfac = 1000.

# ------------------------------------------------------------------------
# Functions

def maybePrint (text):
  if (options.verbose):
    sys.stdout.write(text + '\n')
    sys.stdout.flush()

def getpicksbitmap (infile, dims):
  pickmap = np.zeros(dims)

  with open(infile, 'r') as fp:
    lines = fp.readlines()

  a = lines.index(refloc) + refoffset
  b =  len(lines) - 2

  for i in xrange(a,b):
    sin, srf, pick = lines[i][1:].strip().split('|')[1:-1] 
    sin, srf = int(float(sin.strip()))-1, int(float(srf.strip()))-1
    pick = float(pick.strip())/tfac
    
    if (pick > mintime and pick < maxtime):
      pickmap[sin,srf] = pick

  return np.ma.MaskedArray(pickmap, pickmap == 0)

normfeed = {
	0: {
		'name': 'conventional amp+phase',
		'func': lambda est, obs: est - obs,
	},
	1: {
		'name': 'amplitude-normalized phase only',
		'func': lambda est, obs: np.exp(1j * np.angle(est)) - np.exp(1j * np.angle(obs))
	},
	2: {
		'name': 'amplitude-preserved phase only',
		'func': lambda est, obs: abs(obs)*np.exp(1j * np.angle(est)) - obs
	},
	-1: lambda est, obs: np.log(est/obs).imag,
	-2: lambda est, obs: np.log(est/obs),
	-3: lambda est, obs: np.log(est/obs).real,
}

def calcresid (est, obs, norm, tau = 0., pickmap = None):
  
  diff = normfeed[norm](est, obs)
  if (norm in picknorms and tau):
    taufac = np.exp(pickmap / tau)
    diff *= taufac
  
  return np.nan_to_num(diff * diff.conj()).sum().real

# ------------------------------------------------------------------------

class ModifiedParser(OptionParser):
  def format_epilog(self, formatter):
    return self.epilog

parser = ModifiedParser(usage		= USAGE,
			version		= VERSION,
			description	= DESCRIPTION,
			epilog		= EPILOG)

parser.add_option('-f', '--fazonly', action='store', dest='fazonly',
		help='parameter as in FULLWV; specifies data norm [%default]')

parser.add_option('-p', '--pickfile', action='store', dest='pickfile',
		help='file containing first-arrival picks in ProMax ASCII database output format [%default]')

parser.add_option('-t', '--tau', action='store', dest='tau',
		help='Laplace-domain damping coefficient Tau [%default]')

parser.add_option('-v', '--verbose', action='store_true', dest='verbose',
		help='display additional information')

parser.set_defaults(
			fazonly		= 0,
			pickfile	= 'picks.ascii',
			tau		= 0.,
			verbose		= False,
)

(options, args) = parser.parse_args()

allfiles = glob.glob(gen_obsfilekey)
allfiles.extend(glob.glob(gen_estfilekey))
projnms = np.unique([item.split('.')[0] for item in allfiles])

if (len(args) > 0):
  projnm = args[0]
  if (not projnm in projnms):
    parser.error('No utobs/utest files detected for project \'%s\'!'%(projnm,))
else:
  if (len(projnms) > 1): 
    parser.error(('Please specify a project name; multiple detected:\n' + ['\t\%s\n']*len(projnms))%tuple(projnms))
  elif (len(projnms) == 0):
    parser.error('No utobs/utest files detected!')

  projnm = projnms[0]

norm = int(options.fazonly)
if (not norm in validnorms):
  parser.error('\'%d\' is not a valid norm!'%(norm,))

maybePrint('Processing files for project \'%s\'...'%(projnm,))
lpn = len(projnm)
obsfilekey = proto_obsfilekey%(projnm,)
estfilekey = proto_estfilekey%(projnm,)
obsfilecompose = '%s.utobs'%(projnm,) + '%3.3f'
estfilecompose = '%s.utest'%(projnm,) + '%3.3f'

estfiles = np.array(glob.glob(estfilekey))
sest = np.sort([float(item[lpn+6:]) for item in estfiles])
maybePrint('\tfound %d synthetic data files'%(len(estfiles),))

obsfiles = np.array(glob.glob(obsfilekey))
sobs = np.sort([float(item[lpn+6:]) for item in obsfiles])
maybePrint('\tfound %d observed data files'%(len(obsfiles),))

sortedfreqs = [item for item in sest if item in sobs]
estfiles = [estfilecompose%(item,) for item in sortedfreqs]
obsfiles = [obsfilecompose%(item,) for item in sortedfreqs]
maybePrint('\twith  %d frequencies in common'%(len(sortedfreqs),))

tau = float(options.tau)
if (tau):
  with SEGYFile(obsfiles[0]) as sf:
    dims = (sf.ntr/2, sf.ns)
  if (norm in picknorms):
    if (os.path.isfile(options.pickfile)):
      maybePrint('\treading picks from %s...'%(options.pickfile,))
      pickmap = getpicksbitmap(options.pickfile, dims)
    else:
      parser.error('Laplace-domain (Tau) compensation not implemented for offset! File %s does not exist. Please provide a pick file.'%(options.pickfile,))
  else:
    maybePrint('\tNB: Tau of %3.3f specified; meaningless for fazonly=%d.'%(tau,norm))
else:
  pickmap = None

maybePrint('\tData norm for computation is %d'%(norm,))
maybePrint('\tLaplace time constant is Tau=%3.3f'%(tau,))

residcomposite = 0.
residlist = []
for i in xrange(len(sortedfreqs)):
  sfest = SEGYFile(estfiles[i])
  sfobs = SEGYFile(obsfiles[i])
  utest = sfest[::2] + 1j * sfest[1::2]
  utobs = sfobs[::2] + 1j * sfobs[1::2]

  resid = calcresid(utest, utobs, norm, tau, pickmap)
  residcomposite += resid
  residlist.append(resid)

  print('%g\t%s Hz'%(resid, estfiles[i][12:]))

print('\nTotal residual: %g'%(residcomposite,))

