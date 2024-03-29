all: compile run plot 

compile: edo_2orden.cpp
	g++ -o edo_2orden.exe edo_2orden.cpp

run: edo_2orden.exe
	./edo_2orden.exe

plot: edo_2orden.exe
	python3 plot_edo2orden.py

clean:	
	rm -f  *.dat *.exe 
