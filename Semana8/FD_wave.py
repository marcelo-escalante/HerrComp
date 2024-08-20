# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 14:19:54 2024

@author: Marcelo
"""

import numpy as np
import matplotlib.pyplot as plt

L = 2
c = 300
h0 = 0.1
tf = 0.1

a = 0
b = L
t0 = 0
nx = 100
dx = (b-a)/(nx-1)
dt = 0.3*0.5*dx/c
nt = int((tf-t0)/dt+1)

x = np.linspace(a,b,nx)
t = np.linspace(t0,tf,nt)

u = np.zeros((nx,nt))

# Condiciones de frontera
u[x==a,:] = 0
u[x==b,:] = 0
# u[x==b,:] = h0*np.sin(2*np.pi*100*t) # Función forzada


def f(x):
    return 2*h0/L*x if x <= 0.5*L else -2*h0/L*x+2*h0
f = np.vectorize(f)

# def f(x):
#     return 0.1*x*(L-x)
# f = np.vectorize(f)


# Condiciones iniciales
u[:,0] = f(x)
# u[:,0] = 0 #Para función forzada


cp = dx/dt

#for i in range(1,nx-1):
#    u[i,1] = u[i,0]+0.5*c**2/(cp**2)*(u[i+1,0]+u[i-1,0]-2*u[i,0])

# Usar v(0)=0 para hallar 
u[1:-1,1] = u[1:-1,0]+0.5*c**2/(cp**2)*(u[2:,0]+u[:-2,0]-2*u[1:-1,0])

# Diferencias finitas
for j in range(1,nt-1):
    u[1:-1,j+1] = 2*u[1:-1,j]-u[1:-1,j-1]+c**2/(cp**2)*(u[2:,j]-2*u[1:-1,j]+u[:-2,j])

# print(u[:5,:8])

# for i in range(200):
#     plt.figure()
#     plt.plot(x,u[:,10*i])
#     plt.ylim([-0.11,0.11])
