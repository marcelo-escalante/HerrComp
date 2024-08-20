all: compile run plot 

compile: EscalanteMarcelo_S9C1_EDP.cpp
	g++ -o EscalanteMarcelo_S9C1_EDP.exe EscalanteMarcelo_S9C1_EDP.cpp

run: EscalanteMarcelo_S9C1_EDP.exe
	./EscalanteMarcelo_S9C1_EDP.exe

plot: PLOT_EscalanteMarcelo_S9C1_EDP.exe
	python PLOT_EscalanteMarcelo_S9C1_EDP.py

clean:	
	rm -f  *.dat *.exe 