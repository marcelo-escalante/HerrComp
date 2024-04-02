import numpy as np
import matplotlib.pylab as plt

datos_euler=np.genfromtxt("data1_euler.dat")
datos_leapfrog=np.genfromtxt("data2_leapfrog.dat")
datos_rk4=np.genfromtxt("data3_rk4.dat")


#Euler
x_euler = datos_euler[:,0]
y_euler = datos_euler[:,1]

plt.figure()
plt.plot(x_euler,y_euler)
plt.title("Trayectoria Euler")
plt.grid()
plt.savefig("plot1_euler.pdf")


#Leapfrog
x_leapfrog = datos_leapfrog[:,0]
y_leapfrog = datos_leapfrog[:,1]

plt.figure()
plt.plot(x_leapfrog,y_leapfrog)
plt.title("Trayectoria LeapFrog")
plt.grid()
plt.savefig("plot2_leapfrog.pdf")


#RK4
x_rk4 = datos_rk4[:,0]
y_rk4 = datos_rk4[:,1]

plt.figure()
plt.plot(x_rk4,y_rk4)
plt.title("Trayectoria RK4")
plt.grid()
plt.savefig("plot3_RK4.pdf")




