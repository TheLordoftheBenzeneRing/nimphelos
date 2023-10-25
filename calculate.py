import torch,matplotlib.pyplot as plt, numpy as np,sys
from torch import Tensor

k_J = 1.380649e-23
k_M = 20.83661912e3
c = 2.99792458e8
e = 8.8541878128e-12
h = 6.62607015e-34

def linear(temperature,B,D,exp):

    J = torch.arange(0,100)
    energy = B*J*(J+1) - D*J**2*(J+1)**2
    frequency = torch.diff(energy)
    FWHM = 2*frequency*np.sqrt((2*k_J*temperature*np.log(2))/(1e-27*c**2))

    weight = (2*J + 1) * torch.exp(-energy/(k_M*temperature))
    Q = torch.sum(weight)
    fraction = weight/Q

    nus,lines = [],[]
    mask = fraction[1:]>1e-7
    constant = np.sqrt((1e-27*c**2) / (2*np.pi*k_J*temperature)) * ((16*np.pi**3*(3.33564e-30)**2) / (3*e*h*c**3))
    for (freq,FW,frac) in zip(frequency[mask],FWHM[mask],fraction[1:][mask]):
        nu = torch.linspace(freq-3*FW,freq+3*FW,100)
        line = constant * frac * freq**2 * torch.exp(-(1e-27 * c**2 * (nu-freq)**2) / (2*k_J*temperature*freq**2))
        nus += nu.tolist()
        lines += line.tolist()
    nu2,line2 = [],[]
    mask = fraction[1:len(exp)+1]>1e-7
    for (freq,FW,frac) in zip(exp[mask],FWHM[:len(exp)][mask],fraction[1:len(exp)+1][mask]):
        nu = torch.linspace(freq-3*FW,freq+3*FW,100)
        line = constant * frac * freq**2 * torch.exp(-(1e-27 * c**2 * (nu-freq)**2) / (2*k_J*temperature*freq**2))
        nu2 += nu.tolist()
        line2 += line.tolist()
    plt.plot(nus,np.array(lines)/max(lines))
    plt.plot(nu2,-np.array(line2)/max(line2))

linear(10,57635.96828,0.1835058,torch.tensor([115271,230538,345796,461041,576268,691473,806652,921800,1036912,1151985,1267014,1381995,1496923,1611794,1726603,1841346,1956018,2070616,2185135,2299570]))

def symmetric(temperature,B,D,exp):

    energy = []
    for J in torch.arange(0,100):
        temp = []
        for K in range(-J,J+1):
            if B[1] == B[2]: temp.append(B[1]*J*(J+1) + (B[0]-B[1])*K**2 - D[0]*J**2*(J+1)**2 - D[1]*J*(J+1)*K**2 - D[2]*K**4)
            elif B[0] == B[1]: temp.append(B[1]*J*(J+1) + (B[2]-B[1])*K**2 - D[0]*J**2*(J+1)**2 - D[1]*J*(J+1)*K**2 - D[2]*K**4)
            else: sys.exit("Given molecule is not a symmetric top")
        energy.append(torch.tensor(temp))
    frequency,weight,FWHM = [],[],[]
    for i in range(len(energy)-1):
        for j in range(len(energy[i+1])-2-i):
            temp = energy[i+1][j+1] - energy[i][j]
            frequency.append(temp)
            weight.append((2*(i+1) + 1) * torch.exp(-energy[i+1][j+1]/(k_M*temperature)))
            FWHM.append(2*temp*np.sqrt((2*k_J*temperature*np.log(2))/(1e-27*c**2)))
    frequency,weight,FWHM = torch.tensor(frequency),torch.tensor(weight),torch.tensor(FWHM)
    Q = torch.sum(weight)
    fraction = weight/Q

    nus,lines = [],[]
    constant = np.sqrt((1e-27*c**2) / (2*np.pi*k_J*temperature)) * ((16*np.pi**3*(3.33564e-30)**2) / (3*e*h*c**3))
    mask = fraction>1e-7
    for (freq,FW,frac) in zip(frequency[mask],FWHM[mask],fraction[mask]):
        nu = torch.linspace(freq-3*FW,freq+3*FW,100)
        line = constant * frac * freq**2 * torch.exp(-(1e-27 * c**2 * (nu-freq)**2) / (2*k_J*temperature*freq**2))
        nus += nu.tolist()
        lines += line.tolist()
    nu2,line2 = [],[]
    mask = fraction[:len(exp)]>1e-7
    for (freq,FW,frac) in zip(exp[mask],FWHM[:len(exp)][mask],fraction[:len(exp)][mask]):
        nu = torch.linspace(freq-3*FW,freq+3*FW,100)
        line = constant * frac * freq**2 * torch.exp(-(1e-27 * c**2 * (nu-freq)**2) / (2*k_J*temperature*freq**2))
        nu2 += nu.tolist()
        line2 += line.tolist()
    plt.plot(nus,np.array(lines)/max(lines))
    plt.plot(nu2,-np.array(line2)/max(line2))

symmetric(10,(159140.330,8545.877,8545.877),(0.003,0.163,2.907),torch.tensor([17091.718,34182.727,34183.349,51271.000,51273.980,51274.947,68354.502,68361.035,68364.955,68366.298,85431.224,85442.528,85450.730,85455.622,85457.272,102499.110,102516.573,102530.348,102540.144,102546.024,102547.984,119556.066,119581.168,119601.726,119617.671,119629.100,119635.958,119638.243,136600.150,136634.030,136662.740,136686.190,136704.501,136717.559,136725.396,136728.009,153629.472,153673.424,153711.520,153743.800,153770.224,153790.768,153805.457,153814.273,153817.211]))
plt.show()