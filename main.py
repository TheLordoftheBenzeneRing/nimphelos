import numpy as np,matplotlib.pyplot as plt

def create_table():
    table_element = Element("displayables")
    table_element.select(".display:none").remove_class("display:none")

k_J = 1.380649e-23
k_M = 20.83661912e3
c = 2.99792458e8
e = 8.8541878128e-12
h = 6.62607015e-34

def diatomic():
    temperature,B,D = float(Element("Temp").element.value),float(Element("RotCon1").element.value),float(Element("CenDis1").element.value)

    J = np.arange(0,100)
    energy = B*J*(J+1) - D*J**2*(J+1)**2
    frequency = np.diff(energy)
    FWHM = 2*frequency*np.sqrt((2*k_J*temperature*np.log(2))/(1e-27*c**2))

    # weight = (2*J + 1) * np.exp(-energy/(k_M*temperature))
    # Q = np.sum(weight)
    # fraction = weight/Q

    # nus,lines = [],[]
    # mask = fraction[1:]>1e-7
    # constant = np.sqrt((1e-27*c**2) / (2*np.pi*k_J*temperature)) * ((16*np.pi**3*(3.33564e-30)**2) / (3*e*h*c**3))
    # for (freq,FW,frac) in zip(frequency[mask],FWHM[mask],fraction[1:][mask]):
    #     nu = np.linspace(freq-3*FW,freq+3*FW,100)
    #     line = constant * frac * freq**2 * np.exp(-(1e-27 * c**2 * (nu-freq)**2) / (2*k_J*temperature*freq**2))
    #     nus += nu.tolist()
    #     lines += line.tolist()
    # plt.plot(nus,np.array(lines)/max(lines))
    # display(plt,target="output")
