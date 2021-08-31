#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 13:18:48 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipe, ellipk


def vecPot_loop(I,a,r):
    k = 2.0*np.sqrt(a*r)/(a+r)
    return 1.0e-7*I*(a+r)/r*( (2.0-k*k)*ellipk(k*k) - 2.0*ellipe(k*k) )

def laplaceKernel(r, rp, dtheta):
    return 1.0/np.sqrt(r*r+rp*rp-2*r*rp*np.cos(dtheta))

r = 1.0


rp_range = np.logspace(-8,8,100)

numTheta = 1800
dtheta = 2.0*np.pi/numTheta
theta_range = np.linspace(0.0, 2.0*np.pi - dtheta, numTheta)

#%%
plot_around=[]
for theta in theta_range:
    plot_around.append(laplaceKernel(r, r*(1+0.1), theta+np.pi) * np.cos(theta+np.pi)) # vector potential

plt.figure()
plt.plot(theta_range, plot_around)
plt.xlabel(r"$(\theta + \pi)$ / rad")
plt.ylabel("integrand")
plt.title("vector potential of wire loop:\nLaplace kernel, cos(theta) function")
plt.grid(True)

#%%
blend_a = 4
blend_b = 12
blend_m = np.pi/2.0 # active radius of blending function
def blend(l,a,b, m):
    norm_l = np.abs(l)/m
    if (norm_l < 1.0):
        return (0.5+0.5*np.cos(norm_l**a * np.pi))**b
    else:
        return 0.0

integrand    = []
goodPart     = []
singularPart = []
sumOfTrick   = []
for theta in theta_range:
    integrand.append(   laplaceKernel(r, r*(1+0.1), theta+np.pi) * np.cos(theta+np.pi)                                       )
    goodPart.append(    laplaceKernel(r, r*(1+0.1), theta+np.pi) * np.cos(theta+np.pi)*(1-blend(theta-np.pi, blend_a, blend_b, blend_m)))
    singularPart.append(laplaceKernel(r, r*(1+0.1), theta+np.pi) * np.cos(theta+np.pi)*(  blend(theta-np.pi, blend_a, blend_b, blend_m)))
    sumOfTrick.append(goodPart[-1]+singularPart[-1])

plt.figure()

plt.subplot(3,1,1)
plt.plot(theta_range, integrand, ".", label="original")
plt.plot(theta_range, sumOfTrick, label="good+singular")
plt.legend(loc="upper right")
plt.grid(True)

plt.subplot(3,1,2)
plt.semilogy(theta_range, np.abs(np.divide(np.subtract(integrand,sumOfTrick),integrand)), ".-", label="deviation")
plt.legend(loc="upper right")
plt.grid(True)

plt.subplot(3,1,3)
plt.plot(theta_range, goodPart,     "g", label="good")
plt.plot(theta_range, singularPart, "r", label="singular")
plt.legend(loc="upper right")
plt.grid(True)

plt.tight_layout()

#%%
analytical = []
integral_values = []
for rp in rp_range:
    kernel_eval=[]
    integral = 0.0
    for theta in theta_range:
        kernel_eval.append(laplaceKernel(r, rp, theta) * np.cos(theta)) # vector potential
        integral += kernel_eval[-1]*dtheta
    integral_values.append(integral)
    
    analytical.append(vecPot_loop(1.0, r, rp))
    
    #plt.plot(np.multiply(theta_range, 180.0/np.pi), kernel_eval, label="r'=%.2fm"%(rp))

plt.figure()
plt.loglog(rp_range, np.multiply(1.0e-7, integral_values), ".-", label="quad")
plt.loglog(rp_range, analytical, "-", label="analytical")
plt.xlabel("r' / m")
plt.ylabel("integral value")
plt.title("r=%.2fm"%(r))
plt.legend(loc="upper right")
plt.grid(True)