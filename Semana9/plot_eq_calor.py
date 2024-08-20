# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 17:36:39 2024

@author: LENOVO
"""

import numpy as np
import matplotlib.pyplot as plt

heat0 = np.genfromtxt("heat_0.dat")
heat100 = np.genfromtxt("heat_100.dat")
heat1000 = np.genfromtxt("heat_1000.dat")
heat2500 = np.genfromtxt("heat_2500.dat")

xy = np.genfromtxt("xy.dat")
x = xy[:,0]
y = xy[:,1]


for h in [heat0,heat100,heat1000,heat2500]:
    plt.figure()
    plt.imshow(h.T, cmap='hot', interpolation='nearest',origin="lower",vmin=40,vmax=60,extent=(x.min(), x.max(), y.min(), y.max()), aspect='auto')
    plt.colorbar(label='Temperatura')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Campo escalar T')
    plt.show()
