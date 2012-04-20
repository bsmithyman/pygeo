#!/usr/bin/env python

import numpy as np
from scipy.interpolate import griddata
from pygeo.fastio import readfast

infile = 'vel50'
outfile = 'line15.rayinvr'

y = 0.
olddims = (161,41,61)
newdims = (321,121)
oldori = (477000.,5856.,-2000.)
newori = (477000.,-2000.)
sl = 20
dxold = 100
dxnew = 50
forstr = '%f\t%f\t%f\t%f\n'

oxgrid, ozgrid = np.mgrid[oldori[0]:oldori[0]+olddims[0]*dxold:dxold,
                          oldori[2]:oldori[2]+olddims[2]*dxold:dxold]

nxgrid, nzgrid = np.mgrid[newori[0]:newori[0]+newdims[0]*dxnew:dxnew,
                          newori[1]:newori[1]+newdims[1]*dxnew:dxnew]

fmodel = readfast(infile, olddims)
oldsl = fmodel[:,sl,:]

gridinput = np.column_stack((oxgrid.ravel(), ozgrid.ravel()))
gridoutput = np.column_stack((nxgrid.ravel(), nzgrid.ravel()))

newsl = griddata(gridinput, oldsl.ravel(), gridoutput, method='linear')
newsl.shape = newdims

saveoutput = np.column_stack((nxgrid.T.ravel(), nzgrid.T.ravel(),newsl.T.ravel()))

with open(outfile, 'w') as fp:
  for i in range(saveoutput.shape[0]):
    fp.write(forstr%(saveoutput[i,0], y, -saveoutput[i,1], saveoutput[i,2]))
