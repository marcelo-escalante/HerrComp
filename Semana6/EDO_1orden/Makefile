data = data1_euler.dat data2_rk4.dat data3_analitica.dat data4_error_euler.dat data5_error_rk4.dat
plots = plot1_euler.pdf plot2_RK4.pdf plot3_analitica.pdf plot 4_error.pdf

$(plots): plots_edo1orden.py $(data)
	python plots_edo1orden.py

$(data): edo_1orden.exe
	./edo_1orden.exe

edo_1orden.exe: edo_1orden.cpp
	g++ -o edo_1orden.exe edo_1orden.cpp 