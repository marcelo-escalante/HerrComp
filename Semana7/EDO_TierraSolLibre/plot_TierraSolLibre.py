import numpy as np
import matplotlib.pylab as plt

datos_euler=np.genfromtxt("data1_euler.dat")
datos_leapfrog=np.genfromtxt("data2_leapfrog.dat")
datos_rk4=np.genfromtxt("data3_rk4.dat")


#Euler
xt_euler = datos_euler[:,0]
yt_euler = datos_euler[:,1]
xs_euler = datos_euler[:,2]
ys_euler = datos_euler[:,3]

plt.figure()
plt.plot(xt_euler,yt_euler)
plt.title("Trayectoria Euler Tierra")
plt.grid()
plt.savefig("plot1_eulerTierra.pdf")

plt.figure()
plt.plot(xs_euler,ys_euler)
plt.title("Trayectoria Euler Sol")
plt.grid()
plt.savefig("plot2_eulerSol.pdf")


#Leapfrog
xt_leapfrog = datos_leapfrog[:,0]
yt_leapfrog = datos_leapfrog[:,1]
xs_leapfrog = datos_leapfrog[:,2]
ys_leapfrog = datos_leapfrog[:,3]

plt.figure()
plt.plot(xt_leapfrog,yt_leapfrog)
plt.title("Trayectoria Leapfrog Tierra")
plt.grid()
plt.savefig("plot3_leapfrogTierra.pdf")

plt.figure()
plt.plot(xs_leapfrog,ys_leapfrog)
plt.title("Trayectoria Leapfrog Sol")
plt.grid()
plt.savefig("plot4_leapfrogSol.pdf")


#RK4
xt_rk4 = datos_rk4[:,0]
yt_rk4 = datos_rk4[:,1]
xs_rk4 = datos_rk4[:,2]
ys_rk4 = datos_rk4[:,3]

plt.figure()
plt.plot(xt_rk4,yt_rk4)
plt.title("Trayectoria RK4 Tierra")
plt.grid()
plt.savefig("plot5_RK4Tierra.pdf")

plt.figure()
plt.plot(xs_rk4,ys_rk4)
plt.title("Trayectoria RK4 Sol")
plt.grid()
plt.savefig("plot6_RK4Sol.pdf")




