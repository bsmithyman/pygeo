from scipy import *
from scipy.fftpack import hilbert
import scipy.ndimage as ndimage
import scipy.signal as signal
import scipy.interpolate as interpolate
import pygeo.dipfilt as dipfilt
import pywt
from pylab import *
import sys
gf = ndimage.gaussian_filter
mf = ndimage.median_filter

from pygeo import segyread

globalWavelet = None
globalCounter = 0
globalTraces = 0

def envelope (seismic):
  '''
  Calculates the envelope of the input data (in either the traces or a pseudotrace) and replaces the
  input with the amplitude information.
  '''

  if (seismic.userArg == 'Trace'):
    for i, trace in enumerate(seismic.trData):
      seismic.trData[i,:] = sqrt(trace**2 + hilbert(trace)**2)
  else:
    for i, trace in enumerate(seismic.trHeadersRec[seismic.userArg]):
      seismic.trHeadersRec[seismic.userArg][i,:] = sqrt(trace**2 + hilbert(trace)**2)

def wavedenoise (seismic):
  '''
  Denoises the seismic data (in either the traces or a pseudotrace) based on hard thresholding of the
  representation in the wavelet transform domain.
  '''

  wavelet = 'sym8'
  mode = 'per'
  perc = 0.001
  level = 6
  if (seismic.userArg == 'Trace'):
    for i, trace in enumerate(seismic.trData):
      temp = pywt.wavedec(trace, wavelet, mode=mode, level=level)
      temp = [item * (item > item.max()*perc) for item in temp]
      temp = pywt.waverec(temp, wavelet, mode=mode)
      seismic.trData[i,:] = temp[:len(trace)]
  else:
    for i, trace in enumerate(seismic.trHeadersRec[seismic.userArg]):
      temp = pywt.wavedec(trace, wavelet, mode=mode, level=level)
      temp = [item * (item > item.max()*perc) for item in temp]
      temp = pywt.waverec(temp, wavelet, mode=mode)
      seismic.trHeadersRec[seismic.userArg][i,:] = temp[:len(trace)]

def thresh (seismic):
  '''
  Hard thresholds the seismic data based on a percentage of the maximum amplitude along the trace.
  '''

  perc = float(seismic.userArg)
  for i, trace in enumerate(seismic.trData):
    envelope = sqrt(trace**2 + hilbert(trace)**2)
    seismic.trData[i,:] = (envelope>envelope.max()*perc)*seismic.trData[i,:]

def dipscan (seismic):
  '''
  EXPERIMENTAL
  Calculates the local dip field of the seismic data.
  '''

  [ori,aniso] = dipfilt.orient(float64(seismic.trData),int(seismic.userArg))
  #envelope = seismic.trHeadersRec['ENVELOPE']
  #seismic.trHeadersRec['ORIENTATION'][:] = ori[:]
  #seismic.trHeadersRec['ANISOTROPY'][:] = aniso[:]
  #offsets = seismic.trHeadersRec['OFFSET']
  #offdiff = (offsets - concatenate(([offsets[0]],offsets[:-1])))/10
  #seismic.trHeadersRec['INST_VELOCITY'][:] = (20000 * -tan(ori)) * (aniso > 0.75)
  #goodness = envelope * aniso**2
  #seismic.trHeadersRec['GOODNESS'][:] = goodness
  #seismic.trHeadersRec['FILTERED'][:] = seismic.trData * aniso**2#goodness
  seismic.trData[:] = seismic.trData * aniso**2#goodness

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

def dipemphasis (seismic):
  '''
  EXPERIMENTAL
  Emphasizes local dips in the seismic data by dip-filtering.
  '''

  if seismic.controlState != 2:
    return
  trackprogress(seismic)

  [ori,aniso] = dipfilt.orient(float64(seismic.trData),int(seismic.userArg))
  seismic.trData[:] = seismic.trData * aniso**2

def plotwavecoeffs (seismic):
  '''
  Plots wavelet coefficients for a 1D wavelet decomposition of the seismic data.
  '''

  if seismic.controlState != 2:
    return
  trackprogress(seismic)

  mode = 'per'
  wavelet='sym8'
  #gray()
  jet()

  level = int(seismic.userArg)
  temp = []
  for i in range(level+1):
    temp.append([])

  for i, trace in enumerate(seismic.trData):
    for j,coeffs in enumerate(pywt.wavedec(trace, wavelet, mode=mode, level=level)):
      temp[j].append(sqrt(coeffs**2 + hilbert(coeffs)**2))

  temp = [array(item) for item in temp]

  for i in range(level+1):
    subplot(level+1,1,i+1)
    imshow(temp[i].T,aspect='auto')
    if (i == 0):
      ylabel('Approximation')
      title('Time-domain Wavelet Decomposition (%d levels, %s wavelet)'%(level,wavelet))
    else:
      ylabel('Detail Level %d'%(i,))

  show()

def wavemute (seismic):
  '''
  EXPERIMENTAL
  Mutes wavelet coefficients to improve picking resolution.
  '''

  if seismic.controlState != 2:
    return
  trackprogress(seismic)

  mode = 'sym'
  wavelet = 'sym8'

  level = 6
  mutelevels = [0,1,3,4,5]

  for i, trace in enumerate(seismic.trData):
    coeffs = pywt.wavedec(trace, wavelet, mode=mode, level=level)
    for j in xrange(len(coeffs)):
      try:
        if (mutelevels.index(j)):
          coeffs[j] = zeros(coeffs[j].shape)
      except ValueError:
          pass

    newtrace = pywt.waverec(coeffs, wavelet, mode=mode)
    seismic.trData[i,:] = newtrace[:len(seismic.trData[i])]

def plotdipmodes (seismic):
  '''
  EXPERIMENTAL
  Plots the mode of the orientation field and attempts
  to identify layers by moveout.
  '''

  if seismic.controlState != 2:
    return
  trackprogress(seismic)

  showplot = True

  [ori,aniso] = dipfilt.orient(float64(seismic.trData),float(seismic.userArg))
  dips = []
  for i,trace in enumerate(ori):
    dips.append(trace[argwhere(aniso[i]>0)].mean())

  dips = 180*array(dips)/pi

  if (showplot):
    gray()
    figure()
    subplot(1,3,1)
    imshow(ori.T,aspect='equal')
    subplot(1,3,2)
    imshow(aniso.T,aspect='equal')
    subplot(1,3,3)
    imshow(aniso.T>0.80,aspect='equal')

    figure()
    plot(dips)
    show()
