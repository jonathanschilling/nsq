#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 23:50:33 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import numpy as np
import matplotlib.pyplot as plt

from scipy.special import ellipj, ellipkinc, ellipk, ellipe
from scipy.optimize import root_scalar

eps=0.1

r = 1.0
rp = r+eps

def F(_phi,_ksq):
    return ellipkinc(_phi, _ksq)

#def am(_theta, _ksq):
#    sn,cn,_,_ = ellipj(_theta, _ksq)
#    #return np.sqrt(sn*sn+cn*cn)
#    return np.arctan2(cn,sn)

def theta(_r, _rp, _phi):
    return 2.0/np.sqrt(_r*_r+_rp*_rp-2*_r*_rp) * F(_phi/2.0, -4.0*_r*_rp/(_r*_r+_rp*_rp-2*_r*_rp))

def theta_prime(_r, _rp, _phi):
    return 1.0/np.sqrt(_r*_r+_rp*_rp-2*_r*_rp*(1 - 2.0*np.sin(_phi/2.0)))


def phi(_r, _rp, _theta):
    def diff(_p):
        return _theta-theta(_r, _rp, _p)
    #def diff_prime(_p):
    #    return -theta_prime(_r, _rp, _p)
    #rr=root_scalar(diff, method='newton', bracket=[-np.pi, np.pi], fprime=diff_prime, x0=0.0) # ok, why does this not solve it?
    rr=root_scalar(diff, bracket=[-np.pi, np.pi]) # ok, why does this not solve it?
    return rr.root
    #return 2.0*am(0.5*_theta*np.abs(_r-_rp), -4.0*_r*_rp/(_r*_r+_rp*_rp-2*_r*_rp))

def kernel(_r, _rp, _phi):
    return 1.0/np.sqrt(_r*_r + _rp*_rp - 2.0*_r*_rp*np.cos(_phi))

N = 10+1

plot_leftPhi = -0.2
plot_rightPhi = 0.2

tRange = np.linspace(theta(r, rp, plot_leftPhi), theta(r, rp, plot_rightPhi), N)

tKernels=[]
tPhis = []
tAms = []
for t in tRange:
#    tAms.append(am(t, r, rp))
    p = phi(r, rp, t)
    tPhis.append(p)
    tKernels.append(kernel(r,rp,p))

tKernels = []
for p in tPhis:
    tKernels.append(kernel(r, rp, p))


pRange = np.linspace(plot_leftPhi, plot_rightPhi, N)

pKernels = []
pThetas = []
for p in pRange:
    pKernels.append(kernel(r, rp, p))
    pThetas.append(theta(r, rp, p))

plt.figure(figsize=(8,5))

plt.subplot(2,1,1)
plt.plot(np.multiply(tPhis, 180.0/np.pi), tRange, '.-', label='substitution')
plt.plot(np.multiply(pRange, 180.0/np.pi), pThetas, '.-', label='original')
plt.ylabel(r"$\theta$")
plt.grid(True)
plt.legend(loc='lower right')

plt.subplot(2,1,2)
plt.plot(np.multiply(tPhis, 180.0/np.pi), tKernels, '.-')
plt.plot(np.multiply(pRange, 180.0/np.pi), pKernels, '.-')
plt.ylabel(r"kernel")
plt.grid(True)

plt.xlabel(r"$\varphi$ / deg")

plt.tight_layout()

plt.savefig("nsq_substitution.png")
#%%

# current loop application
# analytical solution
ksq = 4*r*rp/(r*r+rp*rp+2*r*rp)
vecpot_analytical = 1.0e-7 * 4*r/np.sqrt(r*r+rp*rp+2*r*rp)*( (2/ksq-1)*ellipk(ksq) - 2/ksq*ellipe(ksq))
print("analytical: "+str(vecpot_analytical))

# A_x \propto cos(\phi)
def f(_phi):
    return np.cos(_phi)


t0 = theta(r, rp, -np.pi)
t1 = theta(r, rp,  np.pi)

relDevs_p=[]
relDevs_t=[]
allN = np.logspace(1,3,30)
for iN in range(len(allN)):
    N = np.int(allN[iN])
    
    # regular midpoint rule
    dp = 2.0*np.pi/N
    
    pContribs=[]
    for i in range(N):
        p=-np.pi+(i+0.5)*dp
        pContribs.append( f(p) * kernel(r, rp, p) * dp)
    
    vecpot_p = 1.0e-7*np.sum(pContribs)
    relDev_p = (vecpot_p-vecpot_analytical)/vecpot_analytical
#    print("phi grid:   "+str(vecpot_p))
#    print("rel. dev.:  "+str(relDev_p))
    relDevs_p.append(relDev_p)
    
    # transformed midpoint rule
    dt = (t1-t0)/N
    
    ts = []
    tPhis = []
    tContribs = []
    for i in range(np.int(N)):
        t=t0+(i+0.5)*dt
        ts.append(t)
        p = phi(r, rp, t)
        tPhis.append(p)
        #pContribs.append( f(p) * kernel(r, rp, p) * dp)
        tContribs.append( f(p) * dt )
    
    
    vecpot_t = 1.0e-7*np.sum(tContribs)
    relDev_t = (vecpot_t-vecpot_analytical)/vecpot_analytical
#    print("theta grid: "+str(vecpot_t))
#    print("rel. dev.:  "+str(relDev))
    relDevs_t.append(relDev_t)
    
plt.figure()
plt.loglog(allN, np.abs(relDevs_p), label="direct")
plt.loglog(allN, np.abs(relDevs_t), label="transformed")
plt.xlabel("discretization")
plt.ylabel("rel. deviation from analytical solution")
plt.legend(loc='center right')

plt.tight_layout()