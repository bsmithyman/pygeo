#!/usr/bin/env python

# ----------------------------------------------------------------------
from pylab import *
from scipy import *
from pygeo.segyread import *

cdictn = {
       'red':   ((0.0, 1.0, 1.0),
                 (0.10, 1.0, 1.0),
                 (0.30, 1.0, 1.0),
                 (0.40, 0.0, 0.0),
                 (0.60, 0.0, 0.0),
                 (0.65, 0.0, 0.0),
                 (0.75, 1.0, 1.0),
                 (0.95, 1.0, 1.0),
                 (1.0, 1.0, 1.0)),
       'green': ((0.0, 1.0, 1.0),
                 (0.10, 1.0, 1.0),
                 (0.30, 0.0, 0.0),
                 (0.40, 0.0, 0.0),
                 (0.60, 1.0, 1.0),
                 (0.65, 1.0, 1.0),
                 (0.75, 1.0, 1.0),
                 (0.95, 0.0, 0.0),
                 (1.0, 0.0, 0.0)),
       'blue':  ((0.0, 1.0, 1.0),
                 (0.10, 1.0, 1.0),
                 (0.30, 1.0, 1.0),
                 (0.40, 1.0, 1.0),
                 (0.60, 1.0, 1.0),
                 (0.65, 0.0, 0.0),
                 (0.75, 0.0, 0.0),
                 (0.95, 0.0, 0.0),
                 (1.0, 0.33, 0.33)),
}
cdictb = {
       'red':   ((0.0, 1.0, 1.0),
                 (0.5, 1.0, 1.0),
                 (1.0, 0.0, 0.0)),
       'green': ((0.0, 0.0, 0.0),
                 (0.5, 1.0, 1.0),
                 (1.0, 0.0, 0.0)),
       'blue':  ((0.0, 0.0, 0.0),
                 (0.5, 1.0, 1.0),
                 (1.0, 1.0, 1.0))
}
ncomap = matplotlib.colors.LinearSegmentedColormap('seismic', cdictn, 256)
bcomap = matplotlib.colors.LinearSegmentedColormap('seiscomp', cdictb, 256)
nvmi = 1500.
nvma = 6000.
bvm = 500.
asr = 1
aext = (0.,45.,-2.,2.)

# ----------------------------------------------------------------------
sstart = SEGYFile('line10.vp.start')
mstart = sstart.readTraces().T

snew = SEGYFile('line10.vp')
mnew = snew.readTraces().T

sqp = SEGYFile('line10.qp')
mqp = sqp.readTraces().T

# ----------------------------------------------------------------------
figure(1)

ax = subplot(4,1,1)
imshow(mstart, vmin=nvmi, vmax=nvma, cmap=ncomap, aspect=asr, extent=aext)
title('Result from traveltime tomography')
ax.yaxis.set_ticks([-2.0,-1.0,0.0,1.0,2.0])
ylabel('Elevation (km)')
xlabel('Projected line location (km)')
#colorbar()

ax = subplot(4,1,2)
imshow(mnew, vmin=nvmi, vmax=nvma, cmap=ncomap, aspect=asr, extent=aext)
title('Result from waveform tomography')
ax.yaxis.set_ticks([-2.0,-1.0,0.0,1.0,2.0])
ylabel('Elevation (km)')
xlabel('Projected line location (km)')
#colorbar(orientation='horizontal').set_label('Velocity (m/s)')

ax = subplot(4,1,3)
imshow(mnew - mstart, vmin=-bvm, vmax=bvm, cmap=bcomap, aspect=asr, extent=aext)
title('Difference from waveform tomography')
ax.yaxis.set_ticks([-2.0,-1.0,0.0,1.0,2.0])
ylabel('Elevation (km)')
xlabel('Projected line location (km)')
#colorbar(orientation='horizontal').set_label('$\Delta$ Velocity (m/s)')

ax = subplot(4,1,4)
imshow(mqp, vmin=0., vmax=0.15, cmap=cm.gray_r, aspect=asr, extent=aext)
title('Attenuation (1/Q)')
ax.yaxis.set_ticks([-2.0,-1.0,0.0,1.0,2.0])
ylabel('Elevation (km)')
xlabel('Projected line location (km)')
#colorbar(orientation='horizontal')

# ----------------------------------------------------------------------
#ppreA225 = SEGYFile('ppre.start.225.su')
#mppreA225 = array(ppreA225.sNormalize(ppreA225.readTraces(range(1,2363)))).T
#ppreA475 = SEGYFile('ppre.start.475.su')
#mppreA475 = array(ppreA475.readTraces(range(1,2363))).T
#ppreB225 = SEGYFile('ppre.new.225.su')
#mppreB225 = array(ppreB225.readTraces(range(1,2363))).T
#ppreB475 = SEGYFile('ppre.new.475.su')
#mppreB475 = array(ppreB475.readTraces(range(1,2363))).T

#figure(2)

#imshow(mppreA225)

#ax = axes()
#imshow(mnew, vmin=nvmi, vmax=nvma, cmap=ncomap, aspect=asr2, extent=aext)
#ylabel('Elevation (km)')
#xlabel('Projected line location (km)')
#colorbar(orientation='horizontal').set_label('Velocity (m/s)')

show()
