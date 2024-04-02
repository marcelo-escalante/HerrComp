all: compile run plot 

compile: edo_TierraSolFijo.cpp
	g++ -o edo_TierraSolFijo.exe edo_TierraSolFijo.cpp

run: edo_TierraSolFijo.exe
	./edo_TierraSolFijo.exe

plot: edo_TierraSolFijo.exe
	python plot_TierraSolFijo.py

clean:	
	rm -f  *.dat *.exe 