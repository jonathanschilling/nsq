#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 13:18:48 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import numpy as np
import matplotlib.pyplot as plt

def laplaceKernel(r, rp, dtheta):
    return 1.0/np.sqrt(r*r+rp*rp-2*r*rp*np.cos(dtheta))

def laplaceKernel_eps(r, eps, dtheta):
    return 1.0/(r*np.sqrt(eps*eps+2*(1+eps)*(1-np.cos(dtheta))))

r = 1.0


rp_range = np.linspace(0.1, 1.9, 9+1)
theta_p = 180.0 * np.pi/180.0

numTheta = 180
dtheta = 2.0*np.pi/numTheta
theta_range = np.linspace(0.0, 2.0*np.pi - dtheta, numTheta)

plt.figure()
for rp in rp_range:
    
    kernel_eval=[]
    kernel_eps_eval=[]
    integral = 0.0
    for theta in theta_range:
        kernel_eval.append(laplaceKernel(r, rp, theta-theta_p) * np.cos(theta-theta_p)) # vector potential
        eps = (rp-r)/r
        kernel_eps_eval.append(laplaceKernel_eps(r, eps, theta-theta_p) * np.cos(theta-theta_p))
        integral += kernel_eval[-1]
    
    plt.plot(np.multiply(theta_range, 180.0/np.pi), kernel_eval, label="r'=%.2fm"%(rp))
    plt.plot(np.multiply(theta_range, 180.0/np.pi), kernel_eps_eval, ".", label="eps r'=%.2fm"%(rp))
plt.xlabel(r"$\theta$ / deg")
plt.ylabel("Laplace kernel")
plt.title("r=%.2fm theta_prime=%.2f deg"%(r, theta_p*180.0/np.pi))
plt.legend(loc="upper right")
plt.grid(True)

#
##%%
##approach_rp = np.logspace(-3, 3, 100)
#approach_rp = np.linspace(1.0, 100.0, 100)
#singularity=[]
#for rp in approach_rp:
#    #singularity.append(laplaceKernel(1.0, 1.0+rp, 0.0))
#    singularity.append(laplaceKernel(1.0, rp, 0.0))
#
#plt.figure()
#plt.loglog(approach_rp, singularity)
#plt.grid(True)