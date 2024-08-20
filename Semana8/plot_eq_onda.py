# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 23:57:54 2024

@author: LENOVO
"""

import numpy as np
import matplotlib.pylab as plt

data = np.genfromtxt("wave.dat")
x = np.genfromtxt("x.dat")

plt.figure()
plt.ylim([-0.11,0.11])
for i in range(200):
    #plt.figure()
    plt.plot(x,data[:,10*i])
    
plt.savefig("plot_onda.pdf")