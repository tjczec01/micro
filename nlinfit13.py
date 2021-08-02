# -*- coding: utf-8 -*-
"""
Created on Sat Jul 17 09:29:17 2021

@author: tjczec01@gmail.com


"""

from lmfit import minimize, Minimizer, Parameters, Parameter, report_fit
import numpy as np
import scipy.optimize as optimize
import matplotlib.pylab as plt
import math

for i in range(1, 7, 1):
    print("dh{} = \nds{} = \ndg{} = ".format(i, i, i))

def func(P, T, R, dh1, ds1, dg1, dh2, ds2, dg2, dh3, ds3, dg3, dh3t, ds3t, dg3t, dh4, ds4, dg4, dh5, ds5, dg5, dg5t, dh5t, ds5t, dh6, ds6, dg6):
    kB = 1.380662 * 10E23   #  m2 kg/s2 K1
    h = 6.626176 * 10E34  #  m2 kg/s
    Rcal = 1.987  #   cal/K mol
    k1 = math.exp(ds1/R) * math.exp(dh1/(R*T))
    k2 = math.exp(ds2/R) * math.exp(dh2/(R*T))
    k3 = math.exp(ds3/R) * math.exp(dh3/(R*T))
    k4 = math.exp(ds4/R) * math.exp(dh4/(R*T))
    k5 = math.exp(ds5/R) * math.exp(dh5/(R*T))
    k6 = math.exp(ds6/R) * math.exp(dh6/(R*T))
    # Retoh = 



def retoh(P, dh1, ds1, dg1, dh2, ds2, dg2, dh3t, ds3t, dg3t):
    T = 723
    kb = 1.380662 * 10E23   #  m2 kg/s2 K1
    h = 6.626176 * 10E34  #  m2 kg/s
    R = 1.987  #   cal/K mol
    k1 = math.exp(ds1/R) * math.exp((-1*dh1)/(R*T))
    k2 = math.exp(ds2/R) * math.exp((-1*dh2)/(R*T))
    k3 = math.exp(ds3t/R) * math.exp((-1*dh3t)/(R*T))
    
    return ((((kb*T)/h)*(math.exp(ds3t/R) * math.exp((-1*dh3t)/(R*T)))) * math.exp(ds1/R) * math.exp((-1*dh1)/(R*T)) * math.exp(ds2/R) * math.exp((-1*dh2)/(R*T)) * P) / (1 + math.exp(ds1/R) * math.exp((-1*dh1)/(R*T)) * math.exp(ds2/R) * math.exp((-1*dh2)/(R*T)) * P)
    

def rety(P, dh1, ds1, dg1, dh2, ds2, dg2, dh5t, ds5t, dg5t):
    T = 723
    kb = 1.380662 * 10E23   #  m2 kg/s2 K1
    h = 6.626176 * 10E34  #  m2 kg/s
    R = 1.987  #   cal/K mol
    k1 = math.exp(ds1/R) * math.exp((-1*dh1)/(R*T))
    k2 = math.exp(ds2/R) * math.exp((-1*dh2)/(R*T))
    k3 = math.exp(ds5t/R) * math.exp((-1*dh5t)/(R*T))
    
    return (((kb*T)/h)*(math.exp(ds5t/R) * math.exp((-1*dh5t)/(R*T))) * math.exp(ds1/R) * math.exp((-1*dh1)/(R*T)) * P) / (1 + math.exp(ds1/R) * math.exp((-1*dh1)/(R*T)) * math.exp(ds2/R) * math.exp((-1*dh2)/(R*T)) * P)


dh1 = -29.9
ds1 = -26.8 
dg1 = -10.5
dh2 = 5.1
ds2 = 9.9
dg2 = -3.0
dh3 = 30.8
ds3 = 15.9
dg3 = 19.1
dh3t = 51.4
ds3t = 15.9
dg3t = 39.7
dh4 = 9.3
ds4 = 24.5
dg4 = -8.3
dh5 = 34.8
ds5 = 48.5
dg5 = -1.8
dh5t = 70.1
ds5t = 48.5
dg5t = 33.5
dh6 = 9
ds6 = 8.3
dg6 = 1.7

# def func(kd, p0, l0):
#     return 0.5 * (-1 - ((p0 + l0)/kd) + np.sqrt(4 * (l0/kd) + (((l0 - p0)/kd) - 1)**2))

# def residuals(kd, p0, l0, PLP):
#     return PLP - func(kd, p0, l0)

# N=1000
# kd_guess=3.5  # <-- You have to supply a guess for kd
# p0 = np.linspace(0, 10, N)
# l0 = np.linspace(0, 10, N)
# PLP = func(kd_guess, p0, l0) + (np.random.random(N) - 0.5) * 0.1

# kd, cov, infodict, mesg, ier = optimize.leastsq(residuals, kd_guess, args=(p0, l0, PLP), full_output=True)

# print(kd)

# PLP_fit=func(kd, p0, l0)

# plt.plot(p0, PLP, '-b', p0, PLP_fit, '-r')
# plt.show()

# # create data to be fitted
# x = np.linspace(0, 15, 301)
# data = (5.0 * np.sin(2 * x - 0.1) * np.exp(-x * x * 0.025) + np.random.normal(size=len(x), scale=0.2))

# # define objective function: returns the array to be minimized
# def fcn2min(params, x, data):
#     """ model decaying sine wave, subtract data"""
#     amp = params['amp']
#     shift = params['shift']
#     omega = params['omega']
#     decay = params['decay']
#     model = amp * np.sin(x * omega + shift) * np.exp(-x*x*decay)
#     return model - data

# # create a set of Parameters
# params = Parameters()
# params.add('amp',   value= 10,  min=0)
# params.add('decay', value= 0.1)
# params.add('shift', value= 0.0, min=-np.pi/2., max=np.pi/2)
# params.add('omega', value= 3.0)


# # do fit, here with leastsq model
# minner = Minimizer(fcn2min, params, fcn_args=(x, data))
# kws  = {'options': {'maxiter':10}}
# result = minner.minimize()


# # calculate final result
# final = data + result.residual

# # write error report
# report_fit(result)

# # try to plot results
# try:
#     import pylab
#     pylab.plot(x, data, 'k+')
#     pylab.plot(x, final, 'r')
#     pylab.show()
# except:
#     pass