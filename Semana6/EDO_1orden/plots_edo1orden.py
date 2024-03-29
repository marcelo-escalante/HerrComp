# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 10:57:01 2024

@author: Usuario
"""
import numpy as np
import matplotlib.pylab as plt

datos1=np.genfromtxt("y_euler.dat")
datos2=np.genfromtxt("y_rk4.dat")

plt.figure()
plt.plot(datos1[:,0],datos1[:,1])
plt.title("Euler")
plt.grid()
plt.savefig("Euler.pdf")

plt.figure()
plt.plot(datos2[:,0],datos2[:,1])
plt.title("RK4")
plt.grid()
plt.savefig("RK4.pdf")

datos3=np.genfromtxt("error_euler.dat")
datos4=np.genfromtxt("error_rk4.dat")


plt.figure()
plt.plot(datos3[:,0],datos3[:,1], label="Euler")
plt.title("Error")
plt.plot(datos4[:,0],datos4[:,1], label="RK4")
plt.legend(loc = "upper right")
plt.grid()
plt.savefig("Error.pdf")
