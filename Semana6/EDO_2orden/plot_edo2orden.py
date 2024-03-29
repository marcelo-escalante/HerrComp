# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 10:57:01 2024

@author: Usuario
"""

import numpy as np
import matplotlib.pylab as plt


datos1=np.genfromtxt("y_euler.dat")
datos2=np.genfromtxt("y_leap_frog.dat")

datos3=np.genfromtxt("error_euler.dat")
datos4=np.genfromtxt("error_leapfrog.dat")

datos5=np.genfromtxt("y_euler_friccion.dat")
datos6=np.genfromtxt("y_rk4_friccion.dat")

datos7=np.genfromtxt("error_euler_friccion.dat")
datos8=np.genfromtxt("error_rk4_friccion.dat")


plt.figure()
plt.plot(datos1[:,0],datos1[:,1])
plt.title("Euler")
plt.grid()
plt.savefig("Euler.pdf")

plt.figure()
plt.plot(datos2[:,0],datos2[:,1])
plt.title("LeapFrog")
plt.grid()
plt.savefig("LeapFrog.pdf")


plt.figure()
plt.plot(datos3[:,0],datos3[:,1], label="Euler")
plt.plot(datos4[:,0],datos4[:,1], label="LeapFrog")
plt.title("Error")
plt.legend(loc = "upper right")
plt.grid()
plt.savefig("Error.pdf")

plt.figure()
plt.plot(datos5[:,0],datos5[:,1])
plt.title("Euler fricción")
plt.grid()
plt.savefig("Euler_fricción.pdf")

plt.figure()
plt.plot(datos6[:,0],datos6[:,1])
plt.title("RK4 fricción")
plt.grid()
plt.savefig("RK4_fricción.pdf")

plt.figure()
plt.plot(datos7[:,0],datos7[:,1], label="Euler")
plt.title("Error fricción")
plt.plot(datos8[:,0],datos8[:,1], label="RK4")
plt.legend(loc = "upper right")
plt.grid()
plt.savefig("Error_fricción.pdf")
