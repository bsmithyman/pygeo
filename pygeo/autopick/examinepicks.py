#!/usr/bin/env python

from pylab import *
from scipy import *
import scipy.interpolate as interpolate

esmax = 960
infile = 'line5-goldstandard1.pic'

f = open(infile, 'r')
lines = f.readlines()
f.close()

starts = argwhere(array([not line.find('ENSEMBLE NO :') for line in lines]))
metadata = [lines[ind][13:].strip().split() for ind in starts]
metadata = [[int(item[0])]+[int(bit) for bit in item[4:7]] for item in metadata]

shotbitmap = ones((len(starts),esmax))*nan
mx = []
my = []
me = []
offs = []
times = []
for shot in xrange(len(starts)):
  sx = metadata[shot][1]
  sy = metadata[shot][2]
  se = metadata[shot][3]

  ind1 = starts[shot]+1
  try:
    ind2 = starts[shot+1]-1
  except IndexError:
    ind2 = len(lines)
  for rec,time,garbage,rx,ry,re,offset in [line.strip().split() for line in lines[ind1+1:ind2-1]]:
    shotbitmap[shot,int(rec)] = float(time)
    mx.append((float(rx)+float(sx))/20)
    my.append((float(ry)+float(sy))/20)
    me.append((float(re)+float(se))/20)
    offs.append(float(offset)/10)
    times.append(float(time)/10)

#imshow(results)
#show()
