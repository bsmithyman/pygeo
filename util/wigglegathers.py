from pygeo.analysis import *
from pygeo.segyread import SEGYFile

sfs2D = SEGYFile('syn2d-after20-postmod-windowed.ppre.su', isSU=True, endian='Little')
sfs00 = SEGYFile('synth-after00-postmod-windowed.ppre.su', isSU=True, endian='Little')
sfs28 = SEGYFile('synth-after28-postmod-windowed.ppre.su', isSU=True, endian='Little')
sfo = SEGYFile('real-withpicks.ppre.su', isSU=True, endian='Little')

picks = getpicks(sfo)
offsets = calcoffset(sfo)
live = getlive(offsets)
loffs = offsets[live]
lpicks = picks[live]
lrpicks = lpicks/1000. - timereduce(loffs, 4500, 0)

otraces = sfo[:]
s2Dtraces = sfs2D[:]
s00traces = sfs00[:]
s28traces = sfs28[:]

lotraces = otraces[live]
ls2Dtraces = s2Dtraces[live]
ls00traces = s00traces[live]
ls28traces = s28traces[live]

lnotraces = sfo.sNormalize(lotraces)
lns2Dtraces = sfo.sNormalize(ls2Dtraces)
lns00traces = sfo.sNormalize(ls00traces)
lns28traces = sfo.sNormalize(ls28traces)

dto = sfo.trhead[0]['dt']/1000000.
dts2D = sfs2D.trhead[0]['dt']/1000000.
dts00 = sfs00.trhead[0]['dt']/1000000.
dts28 = sfs28.trhead[0]['dt']/1000000.

matplotlib.rcParams.update({'font.size': 20})

fig = plt.figure()
fig.subplots_adjust(hspace=0.7)

ax1 = fig.add_subplot(5,1,1, aspect=100)
wiggle(lnotraces, offsets=loffs, redvel=4500, tshift=0, sampr=dto, color='black', skipt=1, scale=1.0, fill=True, line=False, lwidth=0.1)
#plt.plot(lrpicks, 'r,')
plt.axis([0,960,-0.25,1.25])
plt.ylabel('Time (s)')
plt.xlabel('Channel')
plt.title('Observed field data')
ax1.invert_yaxis()

ax2 = fig.add_subplot(5,1,2, aspect=100)
wiggle(lns28traces, offsets=loffs, redvel=4500, tshift=0, sampr=dts28, color='black', skipt=1, scale=1.0, fill=True, line=False, lwidth=0.0)
#plt.plot(lrpicks, 'r,')
plt.axis([0,960,-0.25,1.25])
plt.ylabel('Time (s)')
plt.xlabel('Channel')
plt.title('2.5D Log-Norm Result')
ax2.invert_yaxis()

ax3 = fig.add_subplot(5,1,3, aspect=100)
wiggle(lnotraces, offsets=loffs, redvel=4500, tshift=0, sampr=dto, color='black', skipt=5, scale=0.25, fill=False, lwidth=1.0)
wiggle(lns00traces, offsets=loffs, redvel=4500, tshift=0, sampr=dts00, color='green', skipt=10, scale=0.25, fill=False, lwidth=1.0)
plt.axis([0,960,-0.25,1.25])
plt.ylabel('Time (s)')
plt.xlabel('Channel')
plt.title('Comparison: Starting model from FAST')
ax3.invert_yaxis()

ax4 = fig.add_subplot(5,1,4, aspect=100)
wiggle(lnotraces, offsets=loffs, redvel=4500, tshift=0, sampr=dto, color='black', skipt=5, scale=0.25, fill=False, lwidth=1.0)
wiggle(lns2Dtraces, offsets=loffs, redvel=4500, tshift=0, sampr=dts2D, color='green', skipt=10, scale=0.25, fill=False, lwidth=1.0)
plt.axis([0,960,-0.25,1.25])
plt.ylabel('Time (s)')
plt.xlabel('Channel')
plt.title('Comparison: 2D Log-Norm Result')
ax4.invert_yaxis()

ax5 = fig.add_subplot(5,1,5, aspect=100)
wiggle(lnotraces, offsets=loffs, redvel=4500, tshift=0, sampr=dto, color='black', skipt=5, scale=0.25, fill=False, lwidth=1.0)
wiggle(lns28traces, offsets=loffs, redvel=4500, tshift=0, sampr=dts28, color='green', skipt=10, scale=0.25, fill=False, lwidth=1.0)
plt.axis([0,960,-0.25,1.25])
plt.ylabel('Time (s)')
plt.xlabel('Channel')
plt.title('Comparison: 2.5D Log-Norm Result')
ax5.invert_yaxis()

plt.show()
