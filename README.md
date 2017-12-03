# RP_Bubble
**A minimal Python script for simulating the Rayleigh-Plesset equation**

I wanted to create own solution data from a Rayleigh-Plesset (RP) equation in order to be able to create plots zooming into the moment of speed reversal at minimal radius. I found [this piece of code](http://pastebin.com/JvzaaUGm) which served me as a starting point.

## What is the Rayleigh-Plesset equation?

The RP equation can describe oscillating acoustic bubbles. A gas bubble in a liquid can oscillate, why is that? Any oscillator is an analogon to a mass-spring system. If you have ever pumped up a bicycle tire you know that air is compressible, and that by compressing air you load potential energy onto a spring. For setting up the equation of motion of a mass on a spring the force equilibrium F=m\*a=k\*x between the spring force and the force arising from inertia has to be written down. In the case of the acoustic bubble the two forces which have to be equilibrated are the pressure difference across the interface (i.e. between inside and outside) versus the force arising from the inertia of the surrounding liquid.

Imagining an ideal gas inside the bubble where the pressure is inversely proportional to the volume and the liquid on the outside with its much higher density and inertia, it is pretty easy to find it plausible that there can be oscillation in the radial motion pattern. If one wants to make the model somewhat realistic, one has to at least account for the following effects:
- Surface tension yields a pressure offset inside the bubble which depends on the surface curvature.
- An ideal gas does not take into account the finite size of atoms or molecules, so one can introduce a Van-der-Waals hard core model.
- The pressure inside the bubble is composed of partial pressures, vapour on the one hand and non-condensible gas on the other.
- Acceleration of the bubble interface generates sound waves emitted into the liquid.

Several different ways of deriving the RP equation are outlined and compared by Leighton [1]

## What's interesting about the Rayleigh-Plesset equation?
Above all, it's the nonlinearity. One can imagine that due to the vapour pressure no real force is needed to enlarge a bubble, but to shrink it beyond the threshold where no more vapour is present means to reduce the space available to the non-condensible gas content, and it is clear that this will let the pressure diverge. Looking at the whole system, one has a large volume of liquid in motion, and when the bubble collapses, the kinetic energy of the heavy liquid in motion is loaded onto the spring of the little amount of gas in the bubble. The result are extraordinary pressure and temerature conditions inside the bubble and immense accelerations in the layers of liquid surrounding it. The adiabatic heating of the bubble content during collapse can even give rise to incandescent plasma, therefore one effect of an acoustically oscillating bubble can be so-called sonoluminescence (SL). In 1993 Löfstedt, Barber & Putterman published an article [2] where they had managed to match simulated RP dynamics with experimental bubble size tracking data.

[1] T. G. Leighton: "Derivation of the Rayleigh-Plesset Equation in Terms of Volume", ISVR Technical Report No. 308, University of Southampton (2007) --> [link](https://eprints.soton.ac.uk/45698/)

[2] Ritva Löfstedt, Bradley P. Barber, Seth Putterman: "Toward a hydrodynamic theory of sonoluminescence", Physics of Fluids A: Fluid Dynamics, volume 5 issue 11 pp. 2911-2928 (1993) --> [link](http://aip.scitation.org/doi/abs/10.1063/1.858700)

