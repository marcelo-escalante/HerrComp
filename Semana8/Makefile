all: compile run plot 

compile: eq_onda.cpp
	g++ -o eq_onda.exe eq_onda.cpp

run: eq_onda.exe
	./eq_onda.exe

plot: plot_eq_onda.exe
	python plot_eq_onda.py

clean:	
	rm -f  *.dat *.exe 
