import torch, matplotlib.pyplot as plt, numpy as np
from torch import Tensor

k_J = 1.38e-23
k_G = 20.83661912
c = 3e8

class molecule():
    def __init__(
            self,
            rot: tuple,  # Rotational constants as a tuple of A,B,C in units of GHz
            dip: tuple,  # Dipole moments along the principle axes in order of A,B,C in units of debye
            mass: float, # Mass of the molecule in amu
            type: int=0  # An integer to designate the type of molecule: 0 -> Diatomic, 1 -> Linear, 2 -> Spherical, 3 -> Symmetric, 4 -> Asymmetric
            ) -> None:
        self.rot,self.dip,self.type = rot,dip,type
        self.mass = 1.66054e-27*mass
    
    def energy(self) -> Tensor:

        if self.type == 0:
            J = torch.arange(0,30)
            return self.rot[1]*J*(J+1)
    
    def get_lines(self,T: float) -> Tensor:
        vs,fs = [],[]

        if self.type == 0:
            J = torch.arange(0,30)
            E = self.energy()
            freq = torch.diff(E)
            P = (2*J + 1) * torch.exp(-E/(k_G*T))
            Q = torch.sum(P)
            frac = P/Q
            FWHM = 2*freq*np.sqrt(2*k_J*T*np.log(2)/(self.mass*c**2))
            for (fr,F,fra) in zip(freq,FWHM,frac[1:]):
                v = torch.linspace(fr-3*F,fr+3*F,101)
                f = fra*np.sqrt((self.mass * c**2) / (2 * np.pi * k_J * T))*torch.exp( - (self.mass * c**2 * (v-fr)**2) / (2 * k_J * T *fr**2) )/fr
                vs += v.tolist()
                fs += f.tolist()
            plt.plot(vs,fs)
            plt.show()

CO = molecule((0,57.897869,0),(0,0.122,0),28.01)
CO.get_lines(100)