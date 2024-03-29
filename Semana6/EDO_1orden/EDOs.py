# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 09:38:38 2024

@author: Usuario
"""

import numpy as np
import matplotlib.pyplot as plt

def euler(f,y0,t):
    y = np.zeros(len(t))
    y[0] = y0
    h = t[1]-t[0]
    for i in range(len(t)-1):
        y[i+1] = y[i]+h*f(t[i],y[i])
    return y

def rk4(f,y0,t):
    y = np.zeros(len(t))
    y[0] = y0
    h = t[1]-t[0]
    for i in range(len(t)-1):
        k1 = h*f(t[i],y[i])
        k2 = h*f(t[i]+h/2,y[i]+k1/2)
        k3 = h*f(t[i]+h/2,y[i]+k2/2)
        k4 = h*f(t[i]+h,y[i]+k3)
        y[i+1] = y[i]+1/6*(k1+2*k2+2*k3+k4)
    return y
        
def f(t,y):
    return -y

t = np.linspace(0,5,10000)

ysol_euler = euler(f,1,t)

ysol_rk4 = rk4(f,1,t)

plt.figure()
plt.plot(t,ysol_euler)
plt.title("Euler")
plt.figure()
plt.title("RK4")
plt.plot(t,ysol_rk4)
        
    