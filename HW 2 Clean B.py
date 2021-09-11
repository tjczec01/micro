# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 06:10:39 2021

@author: tjczec01@gmail.com


"""

import math as mt
import numpy as np
import matplotlib.pylab as plt
import os
import IPython
from IPython.display import display, Latex
from scipy.constants import hbar, Planck, speed_of_light, Avogadro
csc, hsc, NA, hbarsc = speed_of_light, Planck, Avogadro, hbar  # SI units

plt.ioff()
ip = IPython.core.getipython.get_ipython()

clear = os.system('cls')
cwd = os.getcwd()
dir_path = os.path.dirname(os.path.realpath(__file__))
path_fol = r"{}\HomeWork".format(dir_path)

print("\n")
print("Current working directory:\n")
print("{}\n".format(cwd))

try:
    os.mkdir(path_fol)
    print("New Folder was created\n")
    print("Current working directory - Created Folder Path:\n")
    print("{}\n".format(path_fol))
except Exception:
    print("Current working directory - Current Folder Path:\n")
    print("{}\n".format(path_fol))
    pass

font = {'family': 'serif',
        'color': 'black',
        'weight': 'normal',
        'size': 12,
        'ha': 'center',
        'va': 'bottom'}

def k1(a, De):
    t1 = 2.0 * (a**2)
    t2 = De * t1
    return t2

def k2(v, μ):
    pi = 3.141592653589793
    t1 = 4.0 * (pi**2)
    t2 = μ * (v**2)
    return t1 * t2

def ur(m1, m2):
    M1 = m1/NA
    M2 = m2/NA
    MSUM = M1 + M2
    MMul = M1*M2
    MF = MMul / MSUM
    return MF   # kg/molecule

def Vx(k, r, re):
    vv0 = r - re
    vx1 = vv0**2
    vx = 0.5*k*vx1
    return vx

def A(ve, De, μ):
    pi = 3.141592653589793
    a1 = 2.0 * pi * ve 
    a2 = 2.0 * De
    a3 = mt.sqrt(μ/a2)
    af = a1 * a3
    aff = 2*np.pi*ve * np.sqrt(μ/(2*De)) 
    return af  # m^-1

def V(r, re, a, De):
    term1 = r - re
    term2 = term1 * a
    term3 = 1.0 / (np.exp(term2))
    term4 = 1.0 - term3
    term5 = term4**2
    term6 = term5 * De
    return term6

def Vd(delrv, a, De):
    term1 = delrv
    term2 = term1 * a
    term3 = np.exp(-1 * term2)
    term4 = 1.0 - term3
    term5 = term4**2
    term6 = term5 * De
    return term6

def VE(r, re, a, Dev):
    '''The Morse energy for harmonic oscillator at distance r'''
    return Dev * (1 - np.exp(-a*(r-re)))**2

def Vm(r, re, k, D):
    "Morse potential at r given spectroscopic constants"

    alpha = np.sqrt(k/(2*D))
    return D * (1 - np.exp(-alpha*(r - re)))**2

def Evib(n):
    '''Given an (integer) quantum number v, return the vibrational energy in joules.'''
    return hsc*csc * ( (n+0.5)*we - (n+0.5)**2 * wexe )

def En(n, De, h, v):
    
    Xe1 = h * v
    Xe2 = 4.0 * De
    Xe = Xe1/Xe2
    term1 = h * v
    term2 = Xe * term1
    n12 = 1/2
    nterm = n + n12
    Ee = nterm * term1 - (nterm**2) * term2
        
    return Ee

def Eho(ve, n=0):
    """Harmonic oscillator energy (J) of vibrational level v.
    nu is frequency in s^-1"""

    return h*ve*(n+1/2)

def Eaho(ve, n, D0):
    """Anharmonic energy (J) of vibrational level v.
    nu=frequency (s^-1)
    D0=depth (J/molec)"""

    De = D0 + Eho(ve, 0)
    vv = h*ve*(n+1/2) - (h*ve)**2 / (4*De) * (n+1/2)**2 
    return vv  # *(10**9)  # h*nu*(v+1/2) - (h*nu)**2 / (4*De) * (v+1/2)**2

# 1 nm = 10 Å
# 1 J = kg⋅m2⋅s−2
# 8065.73 cm-1 = 96.487 kJ/mol
amutokg = 1.66053892173E-27
we = 2.990e3 * 1e2  # cm^-1 -> m^-1
wexe = 52.8186 * 1e2 # cm^-1 -> m^-1
csc, hsc, NA = speed_of_light, Planck, Avogadro
FAC = 100 * hsc * csc
KFAC = FAC/1000
cmtokj = 96.487 / 8065.73 #  1 cm**-1 = 0.011962587391345855 kJ/mol
kjtocm = 8065.73 / 96.487  #  1 kJ/mol =  83.59395566242084 cm**-1
h = 6.62607004 * (10**-34) # E10-34  # m2 kg / s
Na = 6.0221409 * (10**23)  #  molecules/mol
Hm = 1.00782503207E-3  #  amu
CLm = 34.96885268E-3  #  amu
De = 7.427 * (10**-19) #  J
DekJ = De / 1000 #  kJ
JtoCM = 5.03445 * (10**22)#E22  # cm**-1
DeCm = De*JtoCM  # cm**-1
ve = 9.004 * (10**13)#E13 #  Hz  or # 1/s
re = 1.27455e-8 * 1E-2  # cm
zpe = Eho(ve, 0) 
DF = De + zpe
UR = ur(Hm, CLm)  #  amu
URkg = UR * amutokg # kg
AA = A(ve, De, UR)  #  1/m
AAkg = A(ve, De, UR)
K = 2.0*(DekJ*NA)*(AA**2)
KF = K  #  kJ/mol m**2
# kt = 10**9  #  nm
# ktf = kt**2
# mnm = ktf**-1
# KL = [k1(1.88705191439686622620E10, 7.427E-19), k2(ve, ur(Hm, CLm) * amutokg)]
# K = sum(KL)/len(KL) # k1(1.88705191439686622620 * (10**10), De) * Na * 1E-3  #  kJ/mol*nm^2
rlist = [0]
xstart = 0.08e-9 #  m
xend = 0.6e-9 #  m
xstep = 100 # 1E-13 #  m
delr = [i for i in np.arange(xstart, xend, xstep)] #  m
ogr = [i + re for i in delr]
delrnm = [i * (10**9) for i in np.arange(xstart, xend, xstep)] #  nm
delra = [i * 10 for i in delrnm] #  nm
VRl = [] #  kJ/mol
ens = [i for i in range(0, 8, 1)]
Elist = En(7, De, h, ve)
kk1 = k1(1.88705191439686622620E10, 7.427E-19) 
kk2 = k2(ve, ur(Hm, CLm) * amutokg)

for ii in ens:
    Elist = En(ii, De, h, ve)
for i in delr:
    Vv = Vd(i, AA, De) # * 8.359395566242084E21 * 0.011962587391345855 # 8.359395566242084E21# *83.59395566242084#  * 0.011962587391345855 # * 1000
    VRl.append(Vv/1000)
   
rrb = list(np.linspace(0.08e-9, 0.6e-9, 100000))
rr = list(np.linspace(-0.08e-9, 0.6e-9, 100000))
RVL = [(Na*VE(i, re, AA, DekJ)) for i in rr]
RVLb = [(Na*VE(i, re, AA, DekJ)) for i in rrb]
RXL = [Vx(KF, i, re) for i in rr]
RXLb = [Vx(KF, i, re) for i in rrb]
VRV = np.argmin(np.abs(np.array(VRl) - 45000.0)) # foo(VRl, 399, 401)
XRV = np.argmin(np.abs(np.array(VRl) - 0.1))
Vxl = [] #  kJ/mol
for i in rr:
    Vxv = Vx(KF, i, re) 
    Vxl.append(Vxv)

VRlcm = [(i*NA)*83.59395566242084 for i in VRl]
Vxlcm = [(i*NA)*83.59395566242084 for i in Vxl]
VRlk = [i/1000 for i in VRl]
VRlkcm = [i/1000 for i in VRlcm]
Vxlkcm = [i/1000 for i in Vxlcm]

display(Latex("$\\bf{Question \\ 1}$"))
astring = '$\\bf{A: ' + '{:.5e}'.format(AAkg) + '} \\ \\it{[\\frac{1}{m}]}$'
display(Latex(astring))
# print(kk1, kk2)
# print("K1: {:.3}, K2: {:.3}".format(kk1, k2(ve, ur(Hm, CLm) * amutokg)))
# print("K1: {:.3}, K2: {:.3}".format((k1(1.88705191439686622620 * (10**10), 7.427 * (10**-19))*Na)*mnm  * 0.011962587391345855 * 83.59395566242084 * 1E-3, k2(ur(Hm, CLm) * (10**-10), ve * 83.59395566242084 * 0.011962587391345855 * 1E10)))
# print("A: {:.3}".format(AAkg))
# print("A: {:.20e}".format(AAkg))
# print(cmtokj, kjtocm)
# print(ur(Hm, CLm))
# print(delrnm[VRV], delrnm[-1])
# print(VRl[VRV], Vxl[VRV-1])
# print(Vxl)
# print("URkg: {}".format(URkg))
# print("{:.10f}".format(-1E-10), "{:.8f}".format(-1E-8))
# print(DeCm * 0.011962587391345855, VRl[-1])

DeGraph = DeCm * 0.011962587391345855
xgstart = -0.2
xgend1 = 0.65
xgstep1 = 0.1

xgstart2 = 0.08e-9
xgend2 = 3.5e-10
xgstep2 = 0.5e-10

xgstart3 = -1e-9
xgend3 = 1e-9
xgstep3 = 0.5e-9


figa = plt.figure()
# plt.plot(list(rr), RVL)
plt.plot(rrb, RVLb, color='tab:blue', label=r"$\bf{V(r) = D_{e} * (1 - e^{-a*(r - r_{r})})^2}$")
# plt.plot(rrb, RXLb, color='tab:green', label=r"$\bf{V(r) = k*(r - r_{r})^2}$")
plt.xlabel(r'$\bf{\Delta{R}} \ = \ (r - r_{r})$' + r'$\ _{\it{[nm]}}$')
plt.ylabel(r'$\bf{V(r)}$' +  r'$\ [\frac{kJ}{mol}$]')
plt.legend(loc='best', fontsize='large')
plt.title(label=r"$\bf{Morse} \ \bf{potential} \ (Full)$", fontdict=font)
# plt.xticks(np.arange(xstart, xend, step=xstep))
plt.xlim([0.08e-9, 0.6e-9])
# plt.ylim([0, RVLb[-1]])
plt.grid()
figa.savefig("{}\FullGraphHw2.pdf".format(path_fol), bbox_inches='tight')
figa.savefig("{}\FullGraphHw2.svg".format(path_fol), bbox_inches='tight')
plt.show()
plt.close()

# figb = plt.figure()
# plt.plot(delrnm, VRl, color='tab:blue', label=r"$\bf{V(r) = D_{e} * (1 - e^{-a*(r - r_{r})})^2}$")
# plt.plot(delrnm, Vxl, color='tab:green', label=r"$\bf{V(r) = k*(r - r_{r})^2}$")
# plt.axhline(y=float(DeGraph), color='k', linestyle='--', label=r"$De$")
# plt.xlabel(r'$\bf{\Delta{R}} \ = \ (r - r_{r})$' + r'$\ _{\it{[nm]}}$')
# plt.ylabel(r'$\bf{V(r)}$' +  r'$\ [\frac{kJ}{mol}$]')
# plt.legend(loc='best')
# plt.title(label=r"$\bf{Morse} \ \bf{potential} \ (Full)$", fontdict=font)
# plt.xticks(np.arange(xgstart2, xgend2, step=xgstep2))
# plt.xlim([xgstart2, xgend2])
# plt.ylim(0, VRl[VRV])
# plt.grid()
# figb.savefig("FitGraphHw2.pdf", bbox_inches='tight')
# figb.savefig("FitGraphHw2.svg", bbox_inches='tight')
# plt.show()
# plt.close()

figb = plt.figure()
plt.plot(rr, RVL, color='tab:blue', label=r"$\bf{V(r) = D_{e} * (1 - e^{-a*(r - r_{r})})^2}$")
plt.plot(rr, RXL, color='tab:green', label=r"$\bf{V(r) = k*(r - r_{r})^2}$")
plt.axhline(y=float(DeCm), color='k', linestyle='--', label=r"$De$")
plt.xlabel(r'$Å$')
plt.ylabel(r'$\bf{V(r)}$' +  r'$\ [\frac{1}{cm}$]')
plt.legend(loc='best')
plt.title(label=r"$\bf{Morse} \ \bf{potential} \ (Full)$", fontdict=font)
plt.xticks(np.arange(0.08e-9, 4.25E-10, step=0.5E-10))
plt.xlim([0.08e-9, 4.25E-10])
plt.ylim(0, 500) # VRl[VRV])
plt.grid()
figb.savefig("{}\FitGraphHw2.pdf".format(path_fol), bbox_inches='tight')
figb.savefig("{}\FitGraphHw2.svg".format(path_fol), bbox_inches='tight')
plt.show()
plt.close()

figc = plt.figure()
plt.plot(rr, RVL, color='tab:blue', label=r"$\bf{V(r) = D_{e} * (1 - e^{-a*(r - r_{r})})^2}$")
plt.plot(rr, RXL, color='tab:green', label=r"$\bf{V(r) = k*(r - r_{r})^2}$")
plt.xlabel(r'$Å$')
plt.ylabel(r'$\bf{V(r)}$' +  r'$\ [\frac{1}{cm}$]')
plt.legend(loc='best')
plt.title(label=r"$\bf{Morse} \ \bf{potential}$", fontdict=font)
plt.xticks(np.arange(0.11e-9, 0.15e-9, step=0.01E-9))
plt.xlim([0.11e-9, 0.15e-9])
plt.ylim(0, 20)
plt.grid()
figc.savefig("{}\ZoomedGraphHw2.pdf".format(path_fol), bbox_inches='tight')
figc.savefig("{}\ZoomedGraphHw2.svg".format(path_fol), bbox_inches='tight')
plt.show()
plt.close()