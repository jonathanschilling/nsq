#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 10:10:43 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import numpy as np
import matplotlib.pyplot as plt

eps=0.3
def f(r):
#    return np.divide(1.0, r)
    return np.divide(1.0, np.sqrt(r*r + eps*eps))

def F(a,b):
#    return np.log(b)-np.log(a)
    return np.log(b+np.sqrt(b*b+eps*eps)) - np.log(a+np.sqrt(a*a+eps*eps))

a = 1e-3
b = 1.0

I_analytical = F(a,b)
print("analytical result: %.2e"%I_analytical)

numQuad = 100

# straight-forward trapezoidal quadrature
linAxis=np.linspace(a, b, numQuad)
f_eval = f(linAxis)
I_trapz = 0.0
for i in range(numQuad-1):
    I_trapz += (f_eval[i]+f_eval[i+1])/2.0 * (b-a)/(numQuad-1)
print("\nlinear trapz result: %.2e"%I_trapz)
print("linear trapz rel. dev.: %.2e"%((I_analytical-I_trapz)/I_analytical))

## adapted grid spacing
#adaptAxis = np.linspace(1.0/a, 1.0/b, numQuad)
#ds = (1.0/b - 1.0/a)/(numQuad-1)
#adaptEval = np.divide(f(np.divide(1.0, adaptAxis)), np.multiply(adaptAxis, adaptAxis))
#I_adapt = 0.0
#for i in range(numQuad-1):
#    I_adapt += -1*adaptEval[i] * ds
#    #I_adapt += -1*(adaptEval[i]+adaptEval[i+1])*0.5 * ds
#print("\nadapt trapz result: %.2e"%I_adapt)
#print("adapt trapz rel. dev.: %.2e"%((I_analytical-I_adapt)/I_analytical))

# substitution for 1/r
logAxis=np.linspace(np.log(a), np.log(b), numQuad)
logEval = np.multiply(f(np.exp(logAxis)), np.exp(logAxis))
I_log = 0.0
for i in range(numQuad-1):
    I_log += (logEval[i]+logEval[i+1])/2.0 * (logAxis[i+1]-logAxis[i])
print("\n1/r trapz result: %.2e"%I_log)
print("1/r trapz rel. dev.: %.2e"%((I_analytical-I_log)/I_analytical))

#
#plt.figure()
##plt.loglog(linAxis, f_eval, ".-", label="linear spacing")
##plt.plot(linAxis, f_eval, ".-", label="linear spacing")
#plt.loglog(logAxis, f(logAxis), ".-", label="log spacing")
#plt.xlabel("r")
#plt.ylabel("f")
#plt.grid(True)
#plt.legend(loc="upper right")





