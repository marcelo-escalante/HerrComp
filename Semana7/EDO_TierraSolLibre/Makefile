all: compile run plot 

compile: edo_TierraSolLibre.cpp
	g++ -o edo_TierraSolLibre.exe edo_TierraSolLibre.cpp

run: edo_TierraSolLibre.exe
	./edo_TierraSolLibre.exe

plot: edo_TierraSolLibre.exe
	python plot_TierraSolLibre.py

clean:	
	rm -f  *.dat *.exe 
