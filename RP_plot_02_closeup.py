#!python
"""
using the functions of RP_lib to
- unpickle stored time series results
- add a figure with close-ups at three different zoom levels
  (The close-ups reveal the acceleration taking place at a late time when
   the bubbles are already small and the abrupt deceleration at the turn-around
   point.)

Markus Stokmaier, Weimar, March 2018
"""
import numpy as np
import matplotlib.pyplot as plt

from RP_lib import Bubble, unpickle_thing

twidth1=4e-9
twidth2=0.4e-9
twidth3=0.1e-9
b2=unpickle_thing('pickled_bubble_1.txt')
b3=unpickle_thing('pickled_bubble_2.txt')
bli=[b2,b3]
cli=['b','c']
plt.figure(figsize=[16,4])
ax1=plt.subplot(1,3,1)
for i,b in enumerate(bli):
    idx_rmin=np.argmin(b.R[:,0])
    t_rmin=b.t[idx_rmin]
    idx_t1=np.argmin(np.fabs(b.t-(t_rmin-twidth1/2.)))
    idx_t2=np.argmin(np.fabs(b.t-(t_rmin+twidth1/2.)))
    ax1.plot(1e9*(b.t[idx_t1:idx_t2]-t_rmin), b.R[idx_t1:idx_t2,0]*1e6,cli[i]+'-',lw=2)
ax1.set_ylim(ymin=0)
ax1.set_xlim([-1e9*twidth1/2.,1e9*twidth1/2.])
ax1.set_xticks([-2,-1,0,1,2])
ax1.set_ylabel(r'bubble radius in $\mu$m')
ax1.set_xlabel(r'time window around minimal radius in ns')

ax2=plt.subplot(1,3,2)
for i,b in enumerate(bli):
    idx_rmin=np.argmin(b.R[:,0])
    t_rmin=b.t[idx_rmin]
    idx_t1=np.argmin(np.fabs(b.t-(t_rmin-twidth2/2.)))
    idx_t2=np.argmin(np.fabs(b.t-(t_rmin+twidth2/2.)))
    ax2.plot(1e12*(b.t[idx_t1:idx_t2]-t_rmin), b.R[idx_t1:idx_t2,0]*1e6,cli[i]+'-',lw=2)
ax2.set_ylim(ymin=0)
ax2.set_xlim([-1e12*twidth2/2.,1e12*twidth2/2.])
ax2.set_xticks([-200,-100,0,100,200])
ax2.set_ylabel(r'bubble radius in $\mu$m')
ax2.set_xlabel(r'time window around minimal radius in fs')

ax3=plt.subplot(1,3,3)
for i,b in enumerate(bli):
    idx_rmin=np.argmin(b.R[:,0])
    t_rmin=b.t[idx_rmin]
    idx_t1=np.argmin(np.fabs(b.t-(t_rmin-twidth3/2.)))
    idx_t2=np.argmin(np.fabs(b.t-(t_rmin+twidth3/2.)))
    ax3.plot(1e12*(b.t[idx_t1:idx_t2]-t_rmin), b.R[idx_t1:idx_t2,0]*1e6,cli[i]+'-',lw=2)
ax3.set_ylim(ymin=0)
ax3.set_xlim([-1e12*twidth3/2.,1e12*twidth3/2.])
ax3.set_ylabel(r'bubble radius in $\mu$m')
ax3.set_xlabel(r'time window around minimal radius in fs')

#plt.tight_layout()
#plt.show()
plt.savefig('RPeq_acetone_around_rmin.png',dpi=240,bbox_inches='tight') # see here: http://stackoverflow.com/questions/21288062/second-y-axis-label-getting-cut-off
plt.close()

