#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 17:52:48 2019

@author: Jonathan Schilling (jonathan.schilling@ipp.mpg.de)
"""

import numpy as np
import matplotlib.pyplot as plt

x0 = -1.0
x1 =  1.0
nX =  100+1

y0 = -1.0
y1 =  1.0
nY =  100+1

eps = 0.1

xRange = np.linspace(x0, x1, nX)
yRange = np.linspace(y0, y1, nY)

jet = plt.get_cmap('jet')

def kernel(x,y,eps):
    return np.divide(1.0, np.sqrt(np.add(np.add(np.multiply(x,x), np.multiply(y,y)), np.multiply(eps, eps))))

mg=np.meshgrid(xRange, yRange)

k=kernel(mg[0], mg[1], eps)

plt.figure(1)
plt.clf()
plt.pcolormesh(xRange, yRange, k, cmap=jet)
plt.axis('equal')
plt.xlabel('x / m')
plt.ylabel('y / m')
plt.title('kernel')
plt.colorbar()
plt.tight_layout()

plt.figure(2)
plt.clf()
plt.plot(xRange, kernel(xRange, 0.0, eps), 'k.-')
plt.xlabel('x / m')
plt.title('kernel')
plt.grid(True)
plt.tight_layout()



