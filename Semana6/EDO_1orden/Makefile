all: compile run plot 

compile: edo_1orden.cpp
	g++ -o edo_1orden.exe edo_1orden.cpp

run: edo_1orden.exe
	./edo_1orden.exe

plot: edo_1orden.exe
	python plot_edo1orden.py

clean:	
	rm -f  *.dat *.exe 