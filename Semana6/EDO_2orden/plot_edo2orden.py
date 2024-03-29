import numpy as np
import matplotlib.pylab as plt


datos1=np.genfromtxt("data1_euler.dat")
datos2=np.genfromtxt("data2_leapfrog.dat")
datos3=np.genfromtxt("data3_analitica.dat")
datos4=np.genfromtxt("data4_error_euler.dat")
datos5=np.genfromtxt("data5_error_leapfrog.dat")

datos6=np.genfromtxt("data6_euler_friccion.dat")
datos7=np.genfromtxt("data7_rk4_friccion.dat")
datos8=np.genfromtxt("data8_analitica_friccion.dat")
datos9=np.genfromtxt("data9_error_euler_friccion.dat")
datos10=np.genfromtxt("data10_error_rk4_friccion.dat")


plt.figure()
plt.plot(datos1[:,0],datos1[:,1])
plt.title("Euler")
plt.grid()
plt.savefig("plot1_euler.pdf")

plt.figure()
plt.plot(datos2[:,0],datos2[:,1])
plt.title("LeapFrog")
plt.grid()
plt.savefig("plot2_leapfrog.pdf")

plt.figure()
plt.plot(datos3[:,0],datos3[:,1])
plt.title("Analítica")
plt.grid()
plt.savefig("plot3_analitica.pdf")

plt.figure()
plt.plot(datos4[:,0],datos4[:,1], label="Euler")
plt.plot(datos5[:,0],datos5[:,1], label="LeapFrog")
plt.title("Error")
plt.legend(loc = "upper right")
plt.grid()
plt.savefig("plot4_error.pdf")

plt.figure()
plt.plot(datos6[:,0],datos6[:,1])
plt.title("Euler fricción")
plt.grid()
plt.savefig("plot5_euler_friccion.pdf")

plt.figure()
plt.plot(datos7[:,0],datos7[:,1])
plt.title("RK4 fricción")
plt.grid()
plt.savefig("plot6_RK4_friccion.pdf")

plt.figure()
plt.plot(datos8[:,0],datos8[:,1])
plt.title("Analítica fricción")
plt.grid()
plt.savefig("plot7_analitica_friccion.pdf")

plt.figure()
plt.plot(datos9[:,0],datos9[:,1], label="Euler")
plt.title("Error fricción")
plt.plot(datos10[:,0],datos10[:,1], label="RK4")
plt.legend(loc = "upper right")
plt.grid()
plt.savefig("plot8_error_friccion.pdf")
