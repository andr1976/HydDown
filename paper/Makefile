# Build the IJGGC paper (Elsevier elsarticle). Requires a LaTeX toolchain (latexmk + pdflatex + bibtex).
MAIN = main

.PHONY: all clean
all:
	latexmk -pdf -interaction=nonstopmode $(MAIN).tex

clean:
	latexmk -C $(MAIN).tex
