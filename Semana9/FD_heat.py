# -*- coding: utf-8 -*-
"""
Created on Sun Jul 28 10:30:25 2024

@author: LENOVO
"""

import numpy as np
import matplotlib.pyplot as plt

L = 1
v = 10**-4
T_i = 50
T_ic = 100
T_f = 50
a = 0.2
tf = 2500

dx = 0.01
dt = 0.1*dx**2/v

nx = ny = int(1/dx + 1)
nt = int(tf/dt + 1)

t = np.linspace(0,tf,nt)
x = np.linspace(0,L,nx)
y = np.linspace(0,L,ny)

T = np.zeros((nx,ny,nt))


# Condiciones de frontera

T[x==0,:,:] = T_f
T[x==L,:,:] = T_f
T[:,y==0,:] = T_f
T[:,y==L,:] = T_f


# Condiciones iniciales

T[:,:,0] = T_i

mask_x = (x >= 0.2) & (x <= 0.4)
mask_y = (y >= 0.4) & (y <= 0.6)
T[:,:,0][np.ix_(mask_x, mask_y)] = T_ic
#T[(x>=0.2)*(x<=0.4),(y>=0.2)*(y<=0.4),0] = T_ic

periodica = True
means=[T[:,:,0].mean()]
# Diferencias finitas
for k in range(nt-1):
    T[1:-1,1:-1,k+1] = T[1:-1,1:-1,k] + v*dt/(dx**2)*\
    (T[2:,1:-1,k]+T[:-2,1:-1,k]+T[1:-1,2:,k]+T[1:-1,:-2,k]-4*T[1:-1,1:-1,k])
    
    # if periodica:
    #     T[0,:,k+1] = T[-2,:,k+1]
    #     T[-1,:,k+1] = T[1,:,k+1]
    #     T[:,0,k+1] = T[:,-2,k+1]
    #     T[:,-1,k+1] = T[:,1,k+1]  
        
    means.append(T[:,:,k+1].mean())
    
plt.plot(means)

# for k in range(100):
#     plt.figure()
#     plt.imshow(T[:,:,50*k].T, cmap='hot', interpolation='nearest',origin="lower",vmin=40,vmax=60,extent=(x.min(), x.max(), y.min(), y.max()), aspect='auto')
#     plt.colorbar(label='Temperatura')
#     plt.xlabel('X')
#     plt.ylabel('Y')
#     plt.title('Campo escalar T')
#     plt.show()

ind_100 = int(100/dt)
ind_1000 = int(1000/dt)
ind_2500 = int(2500/dt)
    
for k in [0,ind_100,ind_1000,ind_2500]:
    plt.figure()
    plt.imshow(T[:,:,k].T, cmap='hot', interpolation='nearest',origin="lower",vmin=40,vmax=60,extent=(x.min(), x.max(), y.min(), y.max()), aspect='auto')
    plt.colorbar(label='Temperatura')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Campo escalar T')
    plt.show()
    
    

