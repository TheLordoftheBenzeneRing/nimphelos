import numpy as np,matplotlib.pyplot as plt
from math import ceil
import sys
from js import document

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "cmr10, Computer Modern Serif, DejaVu Serif"
plt.rcParams['font.size'] = 12
plt.rcParams["axes.formatter.use_mathtext"] = True
plt.rcParams["mathtext.fontset"] = "cm"
plt.rcParams['figure.figsize'] = [12, 6]
plt.rcParams["figure.autolayout"] = True

name,energies = [],[]

def create_table(top:str):
    num = 1
    temp = document.getElementById("data_table")
    if top == "diatomic":
        for n in name:
            if "Diatomic" in n: num += 1
        name.append("Diatomic{}".format(str(num)))
        content = temp.innerHTML.replace("</tbody>","<tr>\n<td width=\"5%\"><input type=\"checkbox\" id=\"{}\"></td><td width=\"20%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"15%\">{}</td>\n</tr>\n</tbody>".format("Diatomic{}".format(str(num)),"Diatomic{}".format(str(num)),"",Element("RotCon1").element.value,"",Element("CenDis1").element.value,"","",Element("Temp").element.value))
    elif top == "symmetric":
        for n in name:
            if "Symmetric" in n: num += 1
        name.append("Symmetric{}".format(str(num)))
        content = temp.innerHTML.replace("</tbody>","<tr>\n<td width=\"5%\"><input type=\"checkbox\" id=\"{}\"></td><td class=\"label_column\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"15%\">{}</td>\n</tr>\n</tbody>".format("Symmetric{}".format(str(num)),"Symmetric{}".format(str(num)),Element("RotCon2A").element.value,"",Element("RotCon2C").element.value,Element("CenDis2J").element.value,Element("CenDis2JK").element.value,Element("CenDis2K").element.value,Element("Temp").element.value))
    elif top == "asymmetric":
        for n in name:
            if "Asymmetric" in n: num += 1
        name.append("Asymmetric{}".format(str(num)))
        string = "<tr>\n<td width=\"5%\"><input type=\"checkbox\" id=\"{}\"></td><td class=\"label_column\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"10%\">{}</td>\n<td width=\"15%\">{}</td>\n</tr>\n</tbody>".format("Asymmetric{}".format(str(num)),"Asymmetric{}".format(str(num)),Element("RotCon3A").element.value,Element("RotCon3B").element.value,Element("RotCon3C").element.value,"","","",Element("Temp").element.value)
        if document.querySelector("#dipA").checked: string = string.replace("{}".format(Element("RotCon3A").element.value),"<b>{}</b>".format(Element("RotCon3A").element.value))
        if document.querySelector("#dipB").checked: string = string.replace("{}".format(Element("RotCon3B").element.value),"<b>{}</b>".format(Element("RotCon3B").element.value))
        if document.querySelector("#dipC").checked: string = string.replace("{}".format(Element("RotCon3C").element.value),"<b>{}</b>".format(Element("RotCon3C").element.value))
        content = temp.innerHTML.replace("</tbody>",string)
    else:
        pass
    temp.innerHTML = content

k_J = 1.380649e-23
k_M = 20.83661912e3
c = 2.99792458e8
epsilon = 8.8541878128e-12
h = 6.62607015e-34

def maximum_J(B: float,T: float): return ceil(0.5*np.sqrt(1 - (k_M*T/B)*np.log((1e-30)*k_M*T/B)) - 0.5)

def diatomic():
    create_table("diatomic")

    temperature,B = float(Element("Temp").element.value),float(Element("RotCon1").element.value)
    if Element("CenDis1").element.value == "": D = 0.
    else: D = float(Element("CenDis1").element.value)

    J = np.arange(0,maximum_J(B,temperature)+1)
    energy = B*J*(J+1) - D*J**2*(J+1)**2
    Q = np.sum((2*J + 1) * np.exp(-energy/(k_M*temperature)))
    energies.append((J,energy,Q))

def symmetric():
    create_table("symmetric")

    temperature,A,C,top = float(Element("Temp").element.value),float(Element("RotCon2A").element.value),float(Element("RotCon2C").element.value),str(Element("top").element.value)
    if Element("CenDis2J").element.value == "": D_J = 0.
    else: D_J = float(Element("CenDis2J").element.value)
    if Element("CenDis2JK").element.value == "": D_JK = 0.
    else: D_JK = float(Element("CenDis2JK").element.value)
    if Element("CenDis2K").element.value == "": D_K = 0.
    else: D_K = float(Element("CenDis2K").element.value)

    J,K = [],[]
    for j in range(maximum_J((A+C)/2,temperature)+1):
        for k in range(j+1):
            J.append(j)
            K.append(k)
    J,K = np.array(J),np.array(K)
    if top == "prolate": energy = C*J*(J+1) + (A-C)*K**2 - D_J*J**2*(J+1)**2 - D_JK*J*(J+1)*K**2 - D_K*K**4
    if top == "oblate":  energy = A*J*(J+1) + (C-A)*K**2 - D_J*J**2*(J+1)**2 - D_JK*J*(J+1)*K**2 - D_K*K**4
    Q = np.sum((2*J + 1) * np.exp(-energy/(k_M*temperature)))
    energies.append((J,K,energy,Q))

def asymmetric():
    create_table("asymmetric")

    temperature,A,B,C = float(Element("Temp").element.value),float(Element("RotCon3A").element.value),float(Element("RotCon3B").element.value),float(Element("RotCon3C").element.value)

    J,Ka,Kc,energy = [],[],[],[0.]

    for j in range(maximum_J((A+C)/2,temperature)+1):
        a,c = j,0
        for count in range(2*j+1):
            if count == 0: pass
            elif count%2 == 0: a -= 1
            else:  c += 1
            J.append(j)
            Ka.append(a)
            Kc.append(c)
    J,Ka,Kc = np.array(J),np.array(Ka),np.array(Kc)
    for j in range(1,maximum_J((A+C)/2,temperature)+1):
        H = np.zeros((2*j+1,2*j+1))
        for k in range(2*j+1): H[k][k] = ((A+B)/2)*j*(j+1) + (C - (A+B)/2)*(-k+j)**2
        for k in range(2*j-1):
            K = -j+k
            H[k][k+2] = ((A-B)/4)*np.sqrt((j-K)*(j+K+1)*(j-K-1)*(j+K+2))
            K += 2
            H[k+2][k] = ((A-B)/4)*np.sqrt((j+K)*(j-K+1)*(j+K-1)*(j-K+2))
        for E in np.flip(np.linalg.eigvalsh(H)): energy.append(E)
    energy = np.array(energy)
    Q = np.sum((2*J + 1) * np.exp(-energy/(k_M*temperature)))
    energies.append((J,Ka,Kc,energy,Q))

def plot_spectrum(labels,temperature,designation):
    if len(labels) == 3:
        J,energy,Q = labels
        frequency = np.diff(energy)
        FWHM = 2*frequency*np.sqrt((2*k_J*temperature*np.log(2))/(1e-27*c**2))
        weight = (2*J + 1) * np.exp(-energy/(k_M*temperature))
        fraction = weight/Q
    elif len(labels) == 4:
        J,K,energy,Q = labels
        frequency,FWHM,weight = [],[],[]
        for index,j in enumerate(J):
            if j == 0 or j == K[index]: continue
            freq = energy[index]-energy[J==j-1][K[index]]
            frequency.append(freq)
            FWHM.append(2*freq*np.sqrt((2*k_J*temperature*np.log(2))/(1e-27*c**2)))
            weight.append((2*j + 1) * np.exp(-energy[index]/(k_M*temperature)))
        frequency,FWHM,fraction = np.array(frequency),np.array(FWHM),np.array(weight)/Q
    elif len(labels) == 5:
        J,Ka,Kc,energy,Q = labels
        frequency,FWHM,weight = [],[],[]
        for index,j in enumerate(J):
            if j == 0: continue
            if document.querySelector("#dipA").checked:
                print("A")
                ka,e = Kc[J==j-1][Ka[J==j-1] == Ka[index]],energy[J==j-1][Ka[J==j-1] == Ka[index]]
                if np.size(ka) != 0:
                    for k in range(len(ka)):
                        if ka[k] == Kc[index]+1 or ka[k] == Kc[index]-1:
                            freq = energy[index] - e[k]
                            frequency.append(freq)
                            FWHM.append(2*freq*np.sqrt((2*k_J*temperature*np.log(2))/(1e-27*c**2)))
                            weight.append((2*j + 1) * np.exp(-energy[index]/(k_M*temperature)))
            if document.querySelector("#dipB").checked:
                print("B")
                ka,kc,e = Ka[J==j-1],Kc[J==j-1],energy[J==j-1]
                for k in range(len(ka)):
                    if (ka[k] == Ka[index]+1 or ka[k] == Ka[index]-1) and (kc[k] == Kc[index]+1 or kc[k] == Kc[index]-1):
                        freq = energy[index] - e[k]
                        if freq < 0: continue
                        frequency.append(freq)
                        FWHM.append(2*freq*np.sqrt((2*k_J*temperature*np.log(2))/(1e-27*c**2)))
                        weight.append((2*j + 1) * np.exp(-energy[index]/(k_M*temperature)))
            if document.querySelector("#dipC").checked:
                print("C")
                kc,e = Ka[J==j-1][Kc[J==j-1] == Kc[index]],energy[J==j-1][Kc[J==j-1] == Kc[index]]
                if np.size(kc) != 0:
                    for k in range(len(kc)):
                        if kc[k] == Ka[index]+1 or kc[k] == Ka[index]-1:
                            freq = energy[index] - e[k]
                            frequency.append(freq)
                            FWHM.append(2*freq*np.sqrt((2*k_J*temperature*np.log(2))/(1e-27*c**2)))
                            weight.append((2*j + 1) * np.exp(-energy[index]/(k_M*temperature)))
        frequency,FWHM,fraction = np.array(frequency),np.array(FWHM),np.array(weight)/Q
    else: sys.exit("Not a valid molecule")

    nus,lines = [],[]
    constant = np.sqrt((1e-27*c**2) / (2*np.pi*k_J*temperature)) * ((16*np.pi**3*(3.33564e-30)**2) / (3*epsilon*h*c**3))
    for (freq,FW,frac) in zip(frequency,FWHM,fraction):
        nu = np.linspace(freq-3*FW,freq+3*FW,100)
        line = constant * frac * freq**2 * np.exp(-(1e-27 * c**2 * (nu-freq)**2) / (2*k_J*temperature*freq**2))
        nus += nu.tolist()
        lines += line.tolist()
    nus,lines = np.array(nus),np.array(lines)
    plt.plot(nus/10**6,lines/np.max(lines),label=designation)

def plot():
    temp = document.getElementById("data_table")
    plt.clf()
    for i in range(1,len(name)+1):
        if temp.rows[i].cells[0].querySelector("input[type='checkbox']").checked:
            plot_spectrum(energies[i-1],float(temp.rows[i].cells[-1].innerHTML),name[i-1])
    document.getElementById("output").innerHTML = ""
    plt.xlabel("Frequency (MHz)")
    plt.ylabel("Normalized Arb. Intensity")
    plt.legend(loc="upper right")
    if Element("lbound").element.value != "": plt.xlim(left=float(Element("lbound").element.value))
    if Element("ubound").element.value != "": plt.xlim(right=float(Element("ubound").element.value))
    display(plt,target="output")