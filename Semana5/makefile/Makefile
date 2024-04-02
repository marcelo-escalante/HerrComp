all: resultados.tex

resultados.tex: makedatos.py makedatos1.cpp plotdatos.py plotdatos1.py
	python3 makedatos.py
	g++ makedatos1.cpp -o makedatos1
	./makedatos1
	python3 plotdatos.py
	python3 plotdatos1.py
	pdflatex resultados.tex