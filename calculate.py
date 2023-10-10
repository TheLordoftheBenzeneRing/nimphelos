#!/usr/bin/env python

import torch, matplotlib.pyplot as plt, numpy as np
from torch import Tensor

k_J = 1.380649e-23
k_G = 20.83661912
c = 2.99792458e8
e = 8.8541878128e-12
h = 6.62607015e-34

class molecule():
    def __init__(
            self,
            rotcon: tuple,   # Rotational constants as a tuple of A,B,C in units of GHz
            dipole: tuple,   # Dipole moments along the principle axes in order of A,B,C in units of debye
            mass: float,  # Mass of the molecule in amu
            mtype: int=0,  # An integer to designate the type of molecule: 0 -> Diatomic, 1 -> Linear, 2 -> Spherical, 3 -> Symmetric, 4 -> Asymmetric
            D_J: float=0., # Centrifugal distortion constant
            a: tuple=(0.), # Rotation-Vibration coupling constant
            v: tuple=(0) # vibrational quantum numbers 
            ) -> None:
        self.rot,self.dip,self.type,self.D_J,self.a,self.v = rotcon,dipole,mtype,D_J,a,v
        self.mass = 1.66054e-27*mass
    
    def energy(self) -> Tensor:

        if self.type == 0:
            J = torch.arange(0,30)
            return (self.rot[1] - self.a*(self.v+0.5))*J*(J+1) - self.D_J*J**2*(J+1)**2

    def frequency(self) -> Tensor:
        return torch.diff(self.energy())
    
    def get_lines(self,T: float) -> Tensor:
        vs,fs = [],[]

        if self.type == 0:
            freq = self.frequency()
            P = (2*J + 1) * torch.exp(-E/(k_G*T))
            Q = torch.sum(P)
            frac = P/Q
            FWHM = 2*freq*np.sqrt(2*k_J*T*np.log(2)/(self.mass*c**2))
            constant = ((self.dip[1]*3.33564e-30)**2 * 16 * np.pi**2) * np.sqrt((self.mass * c**2) / (2 * np.pi * k_J * T)) / (3 * e * h * c**3)
            for (fr,F,fra) in zip(freq[frac[1:]>1e-5],FWHM[frac[1:]>1e-5],frac[1:][frac[1:]>1e-5]):
                v = torch.linspace(fr-3*F,fr+3*F,100)
                f = fra*constant*fr**2*torch.exp(-(self.mass*c**2*(v - fr)**2)/(2*k_J*T*fr**2))
                vs += v.tolist()
                fs += f.tolist()
            plt.plot(vs,fs)

CO = molecule(rotcon=(0,57.897869,0),dipole=(0,0.122,0),mass=28.01,mtype=0,D_J=0.)
print(CO.frequency())
CO = molecule(rotcon=(0,57.897869,0),dipole=(0,0.122,0),mass=28.01,mtype=0,D_J=0.000184,a=(0.570800),v=(0))
print(CO.frequency())
plt.show()