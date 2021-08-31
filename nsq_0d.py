#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 14:03:47 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import numpy as np
import matplotlib.pyplot as plt

def kernel(eps, l):
    return 1/np.sqrt(eps*eps+l*l)

def blend(l,a,b, m):
    norm_l = np.abs(l)/m
    if (norm_l < 1.0):
        return (0.5+0.5*np.cos(norm_l**a * np.pi))**b
    else:
        return 0.0

#lRange = np.logspace(-6, 4, 700+1)
lRange = np.linspace(-10, 10, 400+1)

eps = 0.1
blend_a = 1
blend_b = 1
blend_m = 2 # active radius of blending function

integrand    = []
goodPart     = []
singularPart = []
sumOfTrick   = []
for l in lRange:
    integrand.append(   kernel(eps, l)                                       )
    goodPart.append(    kernel(eps,l)*(1-blend(l, blend_a, blend_b, blend_m)))
    singularPart.append(kernel(eps,l)*(  blend(l, blend_a, blend_b, blend_m)))
    sumOfTrick.append(goodPart[-1]+singularPart[-1])

plt.figure()

plt.subplot(2,1,1)
plt.plot(lRange, integrand, ".", label="original")
plt.plot(lRange, sumOfTrick, label="good+singular")
plt.legend(loc="upper right")
plt.grid(True)

plt.subplot(2,1,2)
plt.plot(lRange, goodPart,     "g", label="good")
plt.plot(lRange, singularPart, "r", label="singular")
plt.legend(loc="upper right")
plt.grid(True)
    