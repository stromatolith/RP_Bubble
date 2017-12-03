#!python
"""
starting point: http://pastebin.com/JvzaaUGm

integrating the equation system
Rdot = vdot
vdot = f(R,v)

the Rayleigh-Plesset equation:
    R*ddR + 3dR^2/2 = 1/rho ( ...... )
turns into
    Rdotdot = vdot = -3dR^2/(2R) + 1/(rho*R) * ( ...... )

dynamic viscosity eta = nu * rho   with nu the kinematic viscosity

now producing a nice plot for a bubble in acetone at different driving amplitudes

Markus Stokmaier, Weimar, December 2015
"""
import numpy as np
from numpy import pi, sin
from numpy import array, zeros, linspace
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from cPickle import Pickler, Unpickler

def sphere_volume(r):
    return (4./3.)*pi*r**3

class Bubble(object):
    def __init__(self):
        """
        now own shit
        """
        self.R0 = 0.5e-4       # equilibrium radius
        self.R_start = 0.5e-4       # initial radius
        self.dR_start = 0.          # initial interface speed
        #self.Tamb = 298.       # ambient temperature
        self.p0 = 101325.      # ambient pressure = p_inf, i.e.pressure far away where r--> inf
        #self.pvap = 31.7e2     # vapour pressure in Pa (water) 31.7 hPa at 25 C
        #self.s = 0.0728        # surface tension in N/m (water)
        #self.rho = 998.2071    # density in kg per cubic metre (water)
        #self.nu = 1.5e-3       # viscosity in Pa s (water)  see: https://en.wikipedia.org/wiki/Water_%28data_page%29    https://en.wikipedia.org/wiki/Water_(data_page)
        #self.c = 1435.         # speed of sound in m/s
        self.pvap = 0.093*1e5  # vapour pressure in Pa (acetone)  246 hPa at 20 C
        self.s = (25.2+25*0.1120)*1e-3       # surface tension in N/m (acetone) 25.2mN/m at 25 C with T-coeff -0.1120 according to http://www.surface-tension.de/
        self.rho = 790.        # density in kg per cubic metre (acetone)
        self.nu = 0.4*1e-3        # viscosity in Pa s (acetone) 0.4 mPa s at 10 C according to http://chemister.ru/Database/properties-en.php?dbid=1&id=27
        self.c = 1174.         # speed of sound in m/s
        self.kappa = 1.4 #1.4       # polytropic index 
        self.amp = 0.6*self.p0        # amplitude of acoustic forcing
        self.f = 2e4           # acoustic driving frequency in Hz
        self.t = []  # timegrid
        self.R = []  # will contain the time series of R and dR --> 2D array
        self.Pg = [] # for time series of gas pressure
        self.t_start=0.
        self.dt_coarse=1e-9
        self.dt_fine=1e-12
        self.t_run=19e-6
        self.bubble_radiates=True

    def p_acoustic(self,t):
        return -self.amp*sin(2*pi*self.f*t)

    def dR(self, r, t):
        R = r[0]
        dR = r[1]
        p_gas = (self.p0 + 2*self.s/self.R0 - self.pvap) * (self.R0/R)**(3*self.kappa)
        p_surf = 2*self.s/R
        p_liq = p_gas + self.pvap - p_surf
        p_ext =  self.p0 + self.p_acoustic(t)
        if self.bubble_radiates and (len(self.t)>2):
            Pgdot = (self.Pg[-1]-self.Pg[-2]) / (self.t[-1]-self.t[-2])
            radiation_loss = Pgdot * R/self.c
            ddR = -3*dR**2/(2*R) + 1/(self.rho*R) * ( p_liq - 4*self.nu*dR/R - p_ext + radiation_loss)
        else:
            ddR = -3*dR**2/(2*R) + 1/(self.rho*R) * ( p_liq - 4*self.nu*dR/R - p_ext )
        return np.array([dR, ddR])

    def calculate_Pgas(self,r):
        return (self.p0 + 2*self.s/self.R0 - self.pvap) * (self.R0/r)**(3*self.kappa)

    def integrate_RK4(self):
        func=self.dR
        t=self.t; t.append(self.t_start)
        w=self.R; w.append(array([self.R_start, self.dR_start]))
        self.Pg.append(self.calculate_Pgas(self.R_start))
        while t[-1]<self.t_start+self.t_run:
            if w[-1][0]<0.4*self.R0:
                if w[-1][0]<0.04*self.R0:
                    t.append(t[-1]+0.02*self.dt_fine)
                else:
                    t.append(t[-1]+self.dt_fine)
            else:
                t.append(t[-1]+self.dt_coarse)
            h=t[-1]-t[-2]
            k1=func(w[-1],t[-2])
            k2=func(w[-1]+0.5*h*k1,t[-2]+0.5*h)
            k3=func(w[-1]+0.5*h*k2,t[-2]+0.5*h)
            k4=func(w[-1]+h*k3,t[-1])
            w.append(w[-1]+(h*k1+2*h*k2+2*h*k3+h*k4)/6.)
            self.Pg.append(self.calculate_Pgas(w[-1][0]))

    def list2array(self):
        self.t=array(self.t)
        self.R=array(self.R)
        self.Pg=array(self.Pg)

    def pickle_self(self,suffix):
        ofile=open('pickled_bubble_'+suffix+'.txt','w')
        container=Pickler(ofile)
        container.dump(self)
        ofile.close()


def unpickle_thing(path):
    ifile=open(path,'r')
    container=Unpickler(ifile)
    thing=container.load()
    ifile.close()
    return thing

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

