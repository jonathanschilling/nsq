n#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 15:25:08 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import numpy as np
import matplotlib.pyplot as plt

lRange = np.linspace(-1,1,100+1)

def blend(l,a,b):
    return (0.5+0.5*np.cos(l**a*np.pi))**b

def two_power(l,a,b):
    return (1-np.abs(l)**a)**b

def eta(l):
    return np.exp(-36*np.abs(l**8))




plt.figure()
plt.plot(lRange, blend(lRange, 1,1), "k", label="a=1 b=1")
plt.plot(lRange, blend(lRange, 1,2), label="a=1 b=2")
plt.plot(lRange, blend(lRange, 2,1), label="a=2 b=1")
plt.plot(lRange, blend(lRange, 2,2), label="a=2 b=2")
plt.plot(lRange, blend(lRange, 2,4), label="a=2 b=4")

plt.plot(lRange, eta(lRange), label="$\eta(l)$")
plt.plot(lRange, blend(lRange, 4,12), label="a=4 b=12")
plt.plot(lRange, two_power(lRange, 4,12), label="tp a=4 b=12")


plt.xlabel("l")
plt.ylabel("blending")
plt.title(r"$\left( \frac{1}{2} + \frac{1}{2} \mathrm{cos} \left( l^a \pi \right)  \right)^b$")
plt.grid(True)
plt.legend(loc="upper right")