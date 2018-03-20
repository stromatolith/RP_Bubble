#!python
"""
using the functions of RP_lib to
- simulate a solution to the Rayleigh-Plesset equation
- pickle the time series results

The trajectory comparisons reveal how an increase in driving amplitude enhances
the inertia effect: it leads to bigger bubbles and delays the onset of collapse.

Markus Stokmaier, Weimar, March 2018
"""
import numpy as np
import matplotlib.pyplot as plt

from RP_lib import Bubble

radiation_on=True
if radiation_on:
    suffix=''
else:
    suffix='_norad'

b1=Bubble()
b2=Bubble()
b3=Bubble()
bli=[b1,b2,b3]
amps=[0.85,1.1,1.4]
for i,b in enumerate(bli):
    b.bubble_radiates=radiation_on
    b.amp=amps[i]*b.p0
    b.R0=8e-6
    b.R_start=8e-6
    b.f=18e3
    b.t_run=2./b.f
    b.integrate_RK4()
    b.list2array()
    b.pickle_self(suffix+'_'+str(i))

cli=['k','b','c']
plt.figure(figsize=[8,4])
ax1a=plt.axes()
for i,b in enumerate(bli):
    ax1a.plot(1e6*b.t, b.R[:,0]*1e6,cli[i]+'-',lw=2,label=r'A = {}*$p_0$'.format(amps[i]))
ax1a.set_ylim(ymin=0)
ax1a.set_ylabel(r'bubble radius in $\mu$m')
ax1a.set_xlabel(r'time in $\mu$s')
#plt.legend(loc='upper right')
ax1b=ax1a.twinx()
ax1b.plot(1e6*b.t,b.p_acoustic(b.t)/np.amax(np.fabs(b.p_acoustic(b.t))),'g--',lw=2)
ax1b.set_ylabel('normalised driving pressure signal',color='g')
ax1b.set_ylim([-1,1])
for tl in ax1b.get_yticklabels():
    tl.set_color('g')
#plt.tight_layout()
#plt.show()
plt.savefig('RPeq_acetone'+suffix+'.png',dpi=240,bbox_inches='tight') # see here: http://stackoverflow.com/questions/21288062/second-y-axis-label-getting-cut-off
plt.close()

