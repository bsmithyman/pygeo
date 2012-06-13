
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

from scipy import *
from scipy.fftpack import hilbert
import scipy.ndimage as ndimage
import scipy.signal as signal
import scipy.interpolate as interpolate
import pygeo.dipfilt as dipfilt
import pywt
from pylab import *
import sys
import warnings
gf = ndimage.gaussian_filter
mf = ndimage.median_filter

from pygeo import segyread

globalWavelet = None
globalCounter = 0
globalTraces = 0

method = 'earslope'
refine = 'none'
ltr = 5
showplot = False
derthresh = 5.0
curthresh = 1.0
lmovel = 0#5500
fitness = 100
jumpthresh = 256
earwindow = 40
eardamp = 100
slopegparm = (0,5)
offexclude = 3500
stexclude = 100

def window (ints,wstart=10,function=mean):
  outs = zeros(ints.shape)
  for i in xrange(wstart,len(ints)):
    outs[i] = function(ints[i-wstart:i])

  return outs

def windowstart (ints, function):
  outs = zeros(ints.shape)
  for i in xrange(len(ints)):
    outs[i] = function(ints[:i])
  return outs

def energy (ints):
  return sum(ints**2)

def energyratio (ints, wsmall, eps=1E-5, damp=eardamp):
  b = window(ints,wsmall,energy)
  B = windowstart(ints,energy)
  return b/(B+damp)

def trackprogress (seismic):
  '''
  Prints a counter indicating the number of shot gathers processed, the current SHOTID, and the
  total number of traces processed.
  '''

  global globalCounter
  global globalTraces
  globalCounter += 1
  globalTraces += len(seismic.trData)
  print 'Python Filtering | \tGather: % 5d\tSHOTID: % 5d\tTraces: % 8d\r'%(globalCounter,seismic.trHeadersRec['SHOTID'][0],globalTraces), 
  sys.stdout.flush()

def autopick (seismic):
  '''
  Automatically picks the highest-amplitude seismic signal on each trace, and then
  attempts to regularize the output header information by smoothing with a variety of
  B-spline methods.
  '''

  global globalWavelet
  if seismic.controlState != 2:
    return
  trackprogress(seismic)
  outbuf = seismic.trHeadersRec[seismic.userArg]

  sampr = seismic.reelHeaderRec['SAMP_RATE']/1000.

  if (globalWavelet == None):
    sf = segyread.SEGYFile('wavelet.su', isSU=True)
    globalWavelet = sf.readTraces()
    ns = sf.trhead[0]['ns']
    globalWavelet.shape = (ns,)
    del sf

  function = seismic.trData.copy()

  if (method[:3] == 'env'):
    for i, trace in enumerate(function):
      function[i,:] = sqrt(trace**2 + hilbert(trace)**2)
  elif (method[:3] == 'ear'):
    for i, trace in enumerate(function):
      function[i,:] = energyratio(trace, earwindow)

  #seismic.trHeadersRec['ENVELOPE'][:] = envelope
  #seismic.trHeadersRec['XCOR'][:] = xcor
  #seismic.trHeadersRec['XCOR'][:] = function

  offsets = seismic.trHeadersRec['OFFSET'].copy()

  with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    offsets,monot,unmonot = unique1d(offsets,return_index=True,return_inverse=True)
  function = function[monot]

  pickbuf = offsets.copy()
  envpicks = offsets.copy()
  if (method[3:] == 'max'):
    #for i in range(len(function)):
    #  pickbuf[i] = argmax(function[i,:]) 
    pickbuf = argmax(function,axis=1) 
  elif (method[3:] == 'median'):
    #for i in range(len(function)):
    #  pickbuf[i] = median(flipud(argsort(function[i,:]))[:ltr])
    pickbuf = median(fliplr(argsort(function,axis=1))[:,:ltr],axis=1)
  elif (method[3:] == 'slope'):
    pickbuf = argmax(diff(gf(function,slopegparm),axis=1),axis=1)
  pickbuf[:] = sampr*pickbuf

  if (lmovel):
    pickbuf[:] = pickbuf - lmovel*offsets/10

  #try:
  splc = interpolate.splrep(offsets,pickbuf)
  splder = interpolate.splev(offsets,splc,der=1)
  #except ValueError:
  #  outbuf[:] = 0
  #  print "Failed to produce spline!"
  #  return

  splderenv = sqrt(splder**2 + hilbert(splder)**2)
  fp = argwhere(splderenv < derthresh*splderenv.mean())

  # Curve fitting method
  offsub = offsets[fp].reshape((len(fp),))
  pbsub = pickbuf[fp].reshape((len(fp),))

  splc = interpolate.splrep(offsub,pbsub,k=3,s=fitness*len(offsub))
  envpicks = interpolate.splev(offsets,splc)

  splcur = interpolate.splev(offsub,splc,der=2)
  splcurenv = sqrt(splcur**2 + hilbert(splcur)**2)
  keep = argwhere(splcurenv < curthresh*splcurenv.mean())

  offkep = offsub[keep].reshape((len(keep),))
  pbkep = pbsub[keep].reshape((len(keep),))

  splc = interpolate.splrep(offkep,pbkep,k=1,s=fitness*len(offkep))
  mask = zeros(envpicks.shape)
  mask[fp[keep]] = 1
  envpicks = interpolate.splev(offsets,splc) * mask * (abs(offsets)>offexclude)
  envpicks = envpicks * (envpicks > stexclude)

  # Envelope Autopicking Diagnostics
  if (showplot):
    figure()
    subplot(3,1,1)
    plot(offsets,splder,'b-',label='Function')
    plot(offsets,splderenv,'g-',label='Envelope')
    fill_between(offsets,splderenv,alpha=0.5,facecolor='green')
    plot(offsub,ones(len(offsub))*splderenv.mean()*derthresh,'r-', label='Threshold')
    xlabel('Offset (decimetres)')
    ylabel('Amplitude of First Derivative')
    legend(loc=0)

    subplot(3,1,2)
    plot(offsub,splcur,'b-', label='Function')
    plot(offsub,splcurenv,'g-', label='Envelope')
    fill_between(offsub,splcurenv,alpha=0.5,facecolor='green')
    plot(offsub,ones(len(offsub))*splcurenv.mean()*curthresh,'r-', label='Threshold')
    xlabel('Offset (decimetres)')
    ylabel('Amplitude of Second Derivative')
    legend(loc=0)

    subplot(3,1,3)
    plot(offsets,pickbuf,'r.', label='Initial Picks')
    plot(offsets,envpicks,'gx', label='Filtered Picks')
    xlabel('Offset (decimetres)')
    ylabel('Traveltime (ms)')
    legend(loc=0)
  
  # --------------------------------------------------------------------
  # Pick revision via local cross-correlation of traces
  if (refine == 'xcor'):

    envpicks = envpicks / sampr
    stepcor = []
    corenvelope = []
    shifts = []
    shiftor = len(seismic.trData[0])/2
    for i in xrange(len(seismic.trData)-1):
      stepcor.append(correlate(seismic.trData[i],seismic.trData[i+1],mode=1))
      corenvelope.append(gf(real(sqrt(stepcor[-1]**2 + hilbert(stepcor[-1])**2)),5))
      locs = flipud(corenvelope[-1].argsort())[:15]
      shifts.append(mean(corenvelope[-1][locs]*(shiftor-locs)/sum(corenvelope[-1][locs])))

    stepcor = array(stepcor)
    corenvelope = array(corenvelope)
    shifts = array(shifts)

    if (showplot and False):
      figure()
      subplot(3,1,1)
      imshow(stepcor.T)
      subplot(3,1,2)
      imshow(corenvelope.T)
      subplot(3,1,3)
      plot(shifts)

    # Remove spikes
    shifts = mf(shifts,5)
    newshifts = [0.]
    for i in xrange(len(seismic.trData)-1):
      newshifts.append(shifts[i]+newshifts[-1])
    newshifts = array(newshifts)

    fpicks = envpicks.copy()
    bpicks = envpicks.copy()
    newpicks = envpicks.copy()

    goodpicks = argwhere(envpicks != 0)
    
    #for i in xrange(1,len(envpicks)):
    #  if (((abs(fpicks[i] - fpicks[i-1])>jumpthresh) and (fpicks[i-1] != 0)) or (fpicks[i] == 0)):
    #    fpicks[i] = fpicks[i-1]+shifts[i-1]
    #  l = -(i+1)
    #  if (((abs(bpicks[l] - bpicks[l+1])>jumpthresh) and (bpicks[l+1] != 0)) or (bpicks[l] == 0)):
    #    bpicks[l] = bpicks[l+1]-shifts[l+1]
    #newpicks = (bpicks+fpicks)/2

    # Define some clever function on the forward and backward interpolations
    for i in xrange(len(goodpicks)-1):
      for j in xrange(goodpicks[i]+1,goodpicks[i+1]):
        fpicks[j] = fpicks[j-1] + shifts[j-1]
      for j in xrange(goodpicks[i+1]-1,goodpicks[i],-1):
        bpicks[j] = bpicks[j+1] - shifts[j]
      x0 = goodpicks[i]
      x1 = goodpicks[i+1]
      domain = linspace(0,pi/2,x1-x0+1)
      for j in xrange(goodpicks[i]+1,goodpicks[i+1]):
        newpicks[j] = sqrt((fpicks[j]*cos(domain[j-x0]))**2 + (bpicks[j]*sin(domain[j-x0]))**2)

    # Pick the first pick from the forward prediction or the blended version
    #newpicks = newpicks * (newpicks<fpicks) + fpicks * (fpicks<=newpicks)

    if (showplot):
      figure()
      ax = axes()#subplot(1,4,1)
      gray()
      imshow(seismic.trData.T, aspect='auto')
      axt = ax.axis()
      plot(envpicks, 'g.', label='Initial Pick')
      plot(fpicks,'b-', label='Forward Prediction')
      plot(bpicks,'y-', label='Backward Prediction')
      plot(newpicks,'r-', label='Blended')
      ax.axis(axt)
      ylabel('Sample')
      xlabel('Trace')
      title('Picking output overlaid on seismic gather')

    envpicks[:] = newpicks*sampr

  #ninstphase = []
  #finstphase = []
  #binstphase = []
  #for i,trace in enumerate(seismic.trData):
  #  npick = newpicks[unmonot][i]
  #  fpick = fpicks[unmonot][i]
  #  bpick = bpicks[unmonot][i]
  #  ninstphase.append(arctan2(hilbert(trace)[npick],trace[npick]))
  #  finstphase.append(arctan2(hilbert(trace)[fpick],trace[fpick]))
  #  binstphase.append(arctan2(hilbert(trace)[bpick],trace[bpick]))

  #ninstphase = array(ninstphase)
  #finstphase = array(finstphase)
  #binstphase = array(binstphase)

  if (showplot):
  #  figure()
  #  ax = axes()
  #  plot(offsets[unmonot],180*unwrap(ninstphase)/pi,'r.',label='Blended')
  #  plot(offsets[unmonot],gf(mf(180*unwrap(ninstphase)/pi,5),10),'r-',label='Smooth Blended')
  #  plot(offsets[unmonot],180*unwrap(finstphase)/pi,'b.',label='Forwards')
  #  plot(offsets[unmonot],gf(mf(180*unwrap(finstphase)/pi,5),10),'b-',label='Smooth Forwards')
  #  plot(offsets[unmonot],180*unwrap(binstphase)/pi,'y.',label='Backwards')
  #  plot(offsets[unmonot],gf(mf(180*unwrap(binstphase)/pi,5),10),'y-',label='Smooth Backwards')
  #  legend(loc=0)
  #  ylabel('Unwrapped phase (degrees)')
  #  xlabel('Offset (decimetres)')
  #  title('Instantaneous phase with offset (extracted at pick locations)')

  #  axisrange = axis()
  #  ticks = ax.get_yticks()
  #  ticks = range((int(ticks[0])/360 +1)*360,360,360)
  #  ax.set_yticks(ticks)
  #  grid(True)
  #  axis(axisrange)
    show()

  outbuf[:] = envpicks[unmonot]
  offsets = offsets[unmonot]

  if (lmovel):
    outbuf[:] = outbuf + offsets*lmovel/10

