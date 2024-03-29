import numpy as np
import matplotlib.pylab as plt

datos1=np.genfromtxt("data1_euler.dat")
datos2=np.genfromtxt("data2_rk4.dat")
datos3=np.genfromtxt("data3_analitica.dat")
datos4=np.genfromtxt("data4_error_euler.dat")
datos5=np.genfromtxt("data5_error_rk4.dat")


plt.figure()
plt.plot(datos1[:,0],datos1[:,1])
plt.title("Euler")
plt.grid()
plt.savefig("plot1_euler.pdf")

plt.figure()
plt.plot(datos2[:,0],datos2[:,1])
plt.title("RK4")
plt.grid()
plt.savefig("plot2_RK4.pdf")

plt.figure()
plt.plot(datos3[:,0],datos3[:,1])
plt.title("Analítica")
plt.grid()
plt.savefig("plot3_analitica.pdf")


plt.figure()
plt.plot(datos4[:,0],datos4[:,1], label="Euler")
plt.title("Error")
plt.plot(datos5[:,0],datos5[:,1], label="RK4")
plt.legend(loc = "upper right")
plt.grid()
plt.savefig("plot4_error.pdf")
