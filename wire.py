#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 21:11:23 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import numpy as np
import matplotlib.pyplot as plt


a=-1.0
b=1.0

N = 10

alongX = np.array([np.linspace(a,b,N)]).T
yPos = np.zeros([1,N]).T

eps = 1e-1

L = b-a
Ri = np.sqrt(a*a+eps*eps)
Rf = np.sqrt(b*b+eps*eps)
vecpot_analytical = 2.0e-7*np.arctanh(L/(Ri+Rf))
print("analytical: "+str(vecpot_analytical))



#plt.figure()
#plt.plot(alongX, yPos, '.-')
#plt.scatter(0,eps)
#plt.xlabel("x / m")
#plt.ylabel("y / m")
#plt.tight_layout()
#plt.axis("equal")

def x(i):
    return a+i*(b-a)/(N-1)

def f(x):
    return 1.0/np.sqrt(x*x + eps*eps)


dx = (b-a)/(N-1)

contribs = []
for i in range(N):
    #contribs.append( (f(i)+f(i+1))/2.0 * dx ) # trapezoidal
    contribs.append( f(x(i+0.5)) * dx ) # midpoint
    
vecpot = 1.0e-7*np.sum(contribs)
print("numerical:  "+str(vecpot))
print("rel. dev.:  %.3e"%((vecpot-vecpot_analytical)/vecpot_analytical))


def s(r):
    return np.log(r+np.sqrt(r*r+eps*eps))
def r(s):
    return 0.5*(np.exp(s)-eps*eps*np.exp(-s))

ds = (s(b)-s(a))/N
contribs2 = []
#scale2 = []
#r2 = []
for i in range(N):
    s_i = s(a)+(i+0.5)*ds # midpoint
#    s_i  = s(a)+ i    *ds # trapezoidal
#    s_i1 = s(a)+(i+1)*ds # trapezoidal
    _r = r(s_i)
#    _r1 = r(s_i1)
#    print(f(_r)*np.sqrt(_r*_r+eps*eps))
    contribs2.append( f(_r)*np.sqrt(_r*_r+eps*eps) * ds ) # midpoint
#    contribs2.append( (f(_r)*np.sqrt(_r*_r+eps*eps) + f(_r1)*np.sqrt(_r1*_r1+eps*eps))/2.0 * ds ) # trapezoidal
    #contribs2.append( ds )
#    r2.append(_r)
#    scale2.append(1.0/np.sqrt(_r*_r+eps*eps))

vecpot2 = 1.0e-7*np.sum(contribs2)
print("numerical 2:  "+str(vecpot2))
print("rel. dev. 2:  %.3e"%((vecpot2-vecpot_analytical)/vecpot_analytical))


#plt.plot(r2, scale2, '.-')


