import torch,matplotlib.pyplot as plt, numpy as np
from torch import Tensor

k_J = 1.380649e-23
k_M = 20.83661912e3
c = 2.99792458e8
e = 8.8541878128e-12
h = 6.62607015e-34

def linear(temperature,B):

    J = torch.arange(0,100)
    energy = B*J*(J+1)
    frequency = torch.diff(energy)
    FWHM = 2*frequency*np.sqrt((2*k_J*temperature*np.log(2))/(1e-27*c**2))

    weight = (2*J + 1) * torch.exp(-energy/(k_M*temperature))
    Q = torch.sum(weight)
    fraction = weight/Q

    nus,lines = [],[]
    constant = np.sqrt((1e-27*c**2) / (2*np.pi*k_J*temperature)) * ((16*np.pi**3*(3.33564e-30)**2) / (3*e*h*c**3))
    for (freq,FW,frac) in zip(frequency[fraction[1:]<1e-5],FWHM[fraction[1:]<1e-5],fraction[1:][fraction[1:]>1e-5]):
        nu = torch.linspace(freq-3*FW,freq+3*FW,100)
        line = constant * frac * freq**2 * torch.exp(-(1e-27 * c**2 * (nu-freq)**2) / (2*k_J*temperature*freq**2))
        nus += nu.tolist()
        lines += line.tolist()
    plt.plot(nus,np.array(lines)/max(lines))
    plt.show()

linear(100,57000)