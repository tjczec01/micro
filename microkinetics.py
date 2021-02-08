# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 13:01:57 2021

Github: https://github.com/tjczec01

@author: Travis J Czechorski

E-mail: tjczec01@gmail.com

Website: https://github.com/tjczec01/micro

"""

import math
import numpy as np
import sympy as sp
import scipy as sc
import matplotlib.pyplot as plt
import thermo as tc

"""
Terms list:
    
MASI = most abundant surface intermediate
LHHW = Langmuir−Hinshelwood−Hougen−Watson 
ΔG_i**o = Gibbs free energy changes
K_eq,i = equilibrium constant of the ith elementary reaction
ΔH_i**o = standard enthalpy change of the ith elementary reaction 
ΔS_i**o = entropy change of the ith elementary reaction
DFT = density functional theory
σ_i = stoichiometric numbers of the linear combination of steps that lead to an overall stoichiometric reaction
r_A = rate of adsorption per unit area
σ(T) = the probability that collision of a molecule with the clean surface leads to adsorption (sticking coefficient). Value is between 0 and 1
σ° = the sticking coefficient at a reference temperature T_0
f(θ_r) = a function of surface coverage and takes into account the available surface sites for adsorption
θ_r = reduced coverage and is the ratio of the surface coverage over the surface coverage at surface saturation
m_A = mass of molecule A
v_ij = the stoichiometric coefficient of species j in the elementary step i
v_ij > 0 if species j is a product of step i
v_ij < 0 if species j is a reactant of step i
v_ij = 0 if species j does not participate in step i 
θ_∗ = the concentration of free surface sites
n_j = the number of surface sites occupied by the j-th intermediate
z_i = The reversibility of an elementary step
S_i = Sensitivity of an elementary step
X_RC_i = The degree of rate control for step i
E_i**‡ = enthalpy of formation of the transition state of the rate-determining steps from the gas phase reactants
θ_MASI = the surface coverage of the most abundant surface intermediate
E_app = apparent activation energy
KIE = kinetic isotope effect
r_p = the rate of production of desired product
r_R = the rate of consumption of the reactants
ΔE = heat of reaction
X_s−RC_i = scaled degree of rate control for step i
ε_i is a scaling factor
K_eq_R→TS_i = the equilibrium constant for the formation of the transition state of step i 
r_max_i = The maximum rate, calculated by assuming that the transition state is in equilibrium with the gas phase reactant
β = the reversibility of the overall reaction 
θ_*_b = the fraction of free surface sites 
γ_i = θ_I_i/θ_*
RDS = rate-determining step
ΔBE_TS_DRS = the change in the binding energy of the transition state of the rate-determining step from the reference surface
ΔBE_ref = the change in the binding energy of the descriptor from the reference surface 
χ_TS = the slope in the linear relationship between the transition state and the descriptor

Paper Used: Microkinetic Modeling: A Tool for Rational Catalyst Design

https://doi.org/10.1021/acs.chemrev.0c00394  

"""

main_dict = {}
main_list = []

temps = [450, 425, 400, 375]

tempsk = [i + 273.15 for i in temps]


names = ["diethy_ether", "ethylene", "ethanol", "ethoxythanol", "acetaldehyde", "ethyl_acetate", "acetic_acid", 
         "acetone", "3-hydroxybutanol", "butanone", "1,3-butanediol", "crotonaldehyde", "crotyl_alcohol", "1,3-butadiene", 
         "butyraldehyde", "1-butanol", "1-butene", "(2E)-2-butene", "isobutylene", "carbon_dioxide", "water", "hydrogen", "oxygen"]

cas_nums = ["60-29-7", "74-85-1", "64-17-5", "7518-70-9", "75-07-0", "141-78-6", "64-19-7", 
            "67-64-1", "107-89-1", "78-93-3", "107-88-0", "123-73-9", "6117-91-5", "106-99-0", 
            "123-72-8", "71-36-3", "106-98-9", "107-01-7", "115-11-7", "124-38-9", "7732-18-5", "1333-74-0", "7782-44-7"]

molecular_forms =["C4H10O", "C2H4", "C2H6O", "C4H10O2", "C2H4O", "C4H8O2", "C2H4O2", 
                  "C3H6O", "C4H8O2", "C4H8O", "C4H10O2", "C4H6O", "C4H8O", "C4H6", 
                  "C4H8O", "C4H10O", "C4H8", "C4H8", "C4H8", "CO2", "H2O", "H2", "O2"]

# for i in range(len(names)):
#     main_dict[i] = {'Name': '{}'.format(names[i]), 'CAS': '{}'.format(cas_nums[i]), "Molecular Formula": '{}'.format(molecular_forms[i])}
#     main_list.append({'Name': '{}'.format(names[i]), 'CAS': '{}'.format(cas_nums[i]), "Molecular Formula": '{}'.format(molecular_forms[i])})
    
# print(main_dict)
# print(main_list)

final_dict =    {0: {'Name': 'diethy ether', 'CAS': '60-29-7', 'Molecular Formula': 'C4H10O'}, 
                 1: {'Name': 'ethylene', 'CAS': '74-85-1', 'Molecular Formula': 'C2H4'}, 
                 2: {'Name': 'ethanol', 'CAS': '64-17-5', 'Molecular Formula': 'C2H6O'}, 
                 3: {'Name': 'ethoxythanol', 'CAS': '7518-70-9', 'Molecular Formula': 'C4H10O2'}, 
                 4: {'Name': 'acetaldehyde', 'CAS': '75-07-0', 'Molecular Formula': 'C2H4O'}, 
                 5: {'Name': 'ethyl acetate', 'CAS': '141-78-6', 'Molecular Formula': 'C4H8O2'}, 
                 6: {'Name': 'acetic acid', 'CAS': '64-19-7', 'Molecular Formula': 'C2H4O2'}, 
                 7: {'Name': 'acetone', 'CAS': '67-64-1', 'Molecular Formula': 'C3H6O'}, 
                 8: {'Name': '3-hydroxybutanol', 'CAS': '107-89-1', 'Molecular Formula': 'C4H8O2'}, 
                 9: {'Name': 'butanone', 'CAS': '78-93-3', 'Molecular Formula': 'C4H8O'}, 
                 10: {'Name': '1,3-butanediol', 'CAS': '107-88-0', 'Molecular Formula': 'C4H10O2'}, 
                 11: {'Name': 'crotonaldehyde', 'CAS': '123-73-9', 'Molecular Formula': 'C4H6O'}, 
                 12: {'Name': 'crotyl alcohol', 'CAS': '6117-91-5', 'Molecular Formula': 'C4H8O'}, 
                 13: {'Name': '1,3-butadiene', 'CAS': '106-99-0', 'Molecular Formula': 'C4H6'}, 
                 14: {'Name': 'butyraldehyde', 'CAS': '123-72-8', 'Molecular Formula': 'C4H8O'}, 
                 15: {'Name': '1-butanol', 'CAS': '71-36-3', 'Molecular Formula': 'C4H10O'}, 
                 16: {'Name': '1-butene', 'CAS': '106-98-9', 'Molecular Formula': 'C4H8'}, 
                 17: {'Name': '(2E)-2-butene', 'CAS': '107-01-7', 'Molecular Formula': 'C4H8'}, 
                 18: {'Name': 'isobutylene', 'CAS': '115-11-7', 'Molecular Formula': 'C4H8'}, 
                 19: {'Name': 'carbon dioxide', 'CAS': '124-38-9', 'Molecular Formula': 'CO2'}, 
                 20: {'Name': 'water', 'CAS': '7732-18-5', 'Molecular Formula': 'H2O'}}

final_list = [{'Name': 'diethy ether', 'CAS': '60-29-7', 'Molecular Formula': 'C4H10O'}, 
              {'Name': 'ethylene', 'CAS': '74-85-1', 'Molecular Formula': 'C2H4'}, 
              {'Name': 'ethanol', 'CAS': '64-17-5', 'Molecular Formula': 'C2H6O'}, 
              {'Name': 'ethoxythanol', 'CAS': '7518-70-9', 'Molecular Formula': 'C4H10O2'}, 
              {'Name': 'acetaldehyde', 'CAS': '75-07-0', 'Molecular Formula': 'C2H4O'}, 
              {'Name': 'ethyl acetate', 'CAS': '141-78-6', 'Molecular Formula': 'C4H8O2'}, 
              {'Name': 'acetic acid', 'CAS': '64-19-7', 'Molecular Formula': 'C2H4O2'}, 
              {'Name': 'acetone', 'CAS': '67-64-1', 'Molecular Formula': 'C3H6O'}, 
              {'Name': '3-hydroxybutanol', 'CAS': '107-89-1', 'Molecular Formula': 'C4H8O2'},
              {'Name': 'butanone', 'CAS': '78-93-3', 'Molecular Formula': 'C4H8O'}, 
              {'Name': '1,3-butanediol', 'CAS': '107-88-0', 'Molecular Formula': 'C4H10O2'}, 
              {'Name': 'crotonaldehyde', 'CAS': '123-73-9', 'Molecular Formula': 'C4H6O'}, 
              {'Name': 'crotyl alcohol', 'CAS': '6117-91-5', 'Molecular Formula': 'C4H8O'}, 
              {'Name': '1,3-butadiene', 'CAS': '106-99-0', 'Molecular Formula': 'C4H6'}, 
              {'Name': 'butyraldehyde', 'CAS': '123-72-8', 'Molecular Formula': 'C4H8O'}, 
              {'Name': '1-butanol', 'CAS': '71-36-3', 'Molecular Formula': 'C4H10O'}, 
              {'Name': '1-butene', 'CAS': '106-98-9', 'Molecular Formula': 'C4H8'}, 
              {'Name': '(2E)-2-butene', 'CAS': '107-01-7', 'Molecular Formula': 'C4H8'}, 
              {'Name': 'isobutylene', 'CAS': '115-11-7', 'Molecular Formula': 'C4H8'}, 
              {'Name': 'carbon dioxide', 'CAS': '124-38-9', 'Molecular Formula': 'CO2'}, 
              {'Name': 'water', 'CAS': '7732-18-5', 'Molecular Formula': 'H2O'}]



k_b = [1.380649*(10**-23),           #	J/K Boltzmann constant (m**2 kg)/(s**2 K1)
       8.617333262145*(10**5),      #	eV/K
       1.380649*(10**16)]           #	erg/K

h_p = [6.62607015*(10**-34),         #  J⋅s Plank's constant
       4.135667696*(10**15)]        #  eV⋅s

N_a = 6.02214076*(10**23)           #  mol**1 Avogadro constant

R_g = [8.31446261815324,            #  J/(K⋅mol)  Gas Constant
       8.20573660809596*(10**-5),   #  (m3⋅atm)/(K⋅mol)
       1.987]                       # cal/(K mol)

π = 3.14159265359

def rate_p(k, P, α):
    if len(P) == len(α):
        P_α = [x**y for x, y in zip(P, α)]
        P_α_f = np.prod(P_α)
        r_f = k*P_α_f
        return r_f
    else:
        print("P and α lists must be the same length") 
        
def K_eq(ΔG_i_o, T, R=R_g[0]):
    top = -1.0*ΔG_i_o/(R*T)
    f_val = math.exp(top)
    return f_val

def K_eqb(ΔH_i_o, ΔS_i_o, T, R=R_g[0]):
    v1 = math.exp(ΔS_i_o/R)
    v2 = (-1.0*ΔH_i_o)/(R * T)
    v2b = math.exp(v2)
    f_val = v1*v2b
    return f_val

def S_overall(σ_i, ΔS_i_o):
    if len(σ_i) == len(ΔS_i_o):
        S_ove = [x*y for x, y in zip(σ_i, ΔS_i_o)]
        return sum(S_ove)
    else:
        print("σ_i and ΔS_i_o lists must be the same length") 
        
def H_overall(σ_i, ΔH_i_o):
    if len(σ_i) == len(ΔH_i_o):
        H_ove = [x*y for x, y in zip(σ_i, ΔH_i_o)]
        return sum(H_ove)
    else:
        print("σ_i and ΔH_i_o lists must be the same length") 
        
def K_overall(σ_i, K_eq_i):
    if len(σ_i) == len(K_eq_i):
        K_ove = [x**y for x, y in zip(σ_i, K_eq_i)]
        return np.prod(K_ove)
    else:
        print("σ_i and ΔH_i_o lists must be the same length") 
        
def sticking_coeff(σ_o, T, T_0, E_a, R=R_g[0]):
    v1 = (1/T) - (1/T_0)
    v2 = (-E_a/R)
    v3 = v1*v2
    v4 = math.exp(v3)
    vf = v4*σ_o
    return vf
        
def f_r(θ_r, dissociative_adsorption=False):
    if dissociative_adsorption == True:
        return (1.0 - θ_r)**2 
    elif dissociative_adsorption == False:
        return 1.0 - θ_r

def rate_area(T, θ_r, P, σ_o, E_a, m_a, kb = k_b[0], dissociative_adsorption=False, T_0=298.0):
    v1 = math.sqrt(2.0*π*m_a*kb*T)
    v2 = P/v1
    v3 = f_r(θ_r, dissociative_adsorption)
    v4 = sticking_coeff(σ_o, T, T_0, E_a, R=R_g[0])
    vf = v1*v2*v3*v4
    return vf
    
def k_f(T, θ_r, σ_o, E_a, m_a, kb = k_b[0], dissociative_adsorption=False, T_0=298.0):
    v1 = f_r(θ_r, dissociative_adsorption)
    v2 = sticking_coeff(σ_o, T, T_0, E_a, R=R_g[0])
    v3 = math.sqrt(2.0*π*m_a*kb*T)
    vf = (v1*v2)/v3
    return vf

def r_f(α, ΔH_o_t, ΔS_o_t, T, kb=k_b[0], H_p=h_p[0], R=R_g[0]):
    v1 = (kb*T)/H_p
    v2 = math.exp(ΔS_o_t/R)
    v3a = (-1.0*ΔH_o_t)/(R * T)
    v3 = math.exp(v3a)
    v4 = np.prod(α)
    vf = v1*v2*v3*v4
    return vf
    
def k_f_2(ΔH_o_t, ΔS_o_t, T, kb=k_b[0], H_p=h_p[0], R=R_g[0]):
    v1 = (kb*T)/H_p
    v2 = math.exp(ΔS_o_t/R)
    v3a = (-1.0*ΔH_o_t)/(R * T)
    v3 = math.exp(v3a)
    vf = v1*v2*v3
    return vf

# Both the forward and reverse rate constant can be estimated by transition state theory if the values for ΔH°, ΔS°, ΔH°‡, and ΔS°‡ are available

def rate_i(k_i_f, k_i_b, P_j, v_i_j, θ_j):
    v1a = [x**(-y) for x, y in zip(P_j, v_i_j)]
    v1b = [x**(-y) for x, y in zip(θ_j, v_i_j)]
    v1 = k_i_f*np.prod(v1a)*np.prod(v1b)
    
    v2a = [x**y for x, y in zip(P_j, v_i_j)]
    v2b = [x**y for x, y in zip(θ_j, v_i_j)]
    v2 = k_i_b*np.prod(v2a)*np.prod(v2b)
    
    vf = v1 - v2
    return vf

def z_tot(z_i, σ_i):
    v1 = [x**y for x,y in zip(z_i, σ_i)]
    return np.prod(v1)

def E_app(X_RC, E_t, n, θ_MASI, H_MASI_g):
    # n is the number of surface sites involved in the rate-determining step
    # H_MASI g is the enthalpy of formation of MASI from the gas phase reactants and products 
    # θ_MASI is the surface coverage of the most abundant surface intermediate
    v1a = [x*y for x, y in zip(X_RC, E_t)]
    v1 = sum(v1a)
    v2 = n*θ_MASI*H_MASI_g
    vf = v1 - v2
    return vf

def KIE(G_i_H_o, G_i_D, X_i, T, R=R_g[0]):
    v1a = [(z*(x - y))/(R * T) for x, y, z in zip(G_i_H_o, G_i_D, X_i)]
    v1b = [math.exp(x_i) for x_i in v1a]
    vf = np.prod(v1b)
    return vf

def r_overall(r_max_i, X_RC_i, β, θ_o=1, n=0):
    # β = z1*z2*z3*z4
    # θ_o is the fraction of free surface sites
    r_ove = r_max_i*X_RC_i*(1.0 - β)*(θ_o**n)
    return r_ove

def r_s_RDS(ν_k_MASI, T, P_k, r_max_RDS, K_ads_MASI, ΔBE_TS_DRS, ΔBE_MASI, n, R=R_g[0]):
    v1a = (-1.0*ΔBE_TS_DRS)/(R * T)
    v1b = r_max_RDS*math.exp(v1a)
    v2a = (-1.0*ΔBE_MASI)/(R * T)
    v2b = math.exp(v2a)
    v2c = 1.0 + K_ads_MASI*np.prod([x**(-y) for x, y in zip(P_k, ν_k_MASI)])*v2b
    v3 = v1b/(v2c**n)
    return v3

def r_s_RDSb(ν_k_MASI, T, P_k, r_max_RDS, K_ads_MASI, ΔBE_ref, n, χ_TS, χ_MASI, R=R_g[0]):
    v1a = (-1.0*χ_TS*ΔBE_ref)/(R * T)
    v1b = r_max_RDS*math.exp(v1a)
    v2a = (-1.0*χ_MASI*ΔBE_ref)/(R * T)
    v2b = math.exp(v2a)
    v2c = 1.0 + K_ads_MASI*np.prod([x**(-y) for x, y in zip(P_k, ν_k_MASI)])*v2b
    v3 = v1b/(v2c**n)
    return v3
