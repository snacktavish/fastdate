LATEX = latex
DVIPDF = dvipdf
BIBTEX = bibtex
PDFLATEX = pdflatex

all: dating

dating: dating.tex
	$(LATEX) $+
	$(BIBTEX) $@
	$(DVIPDF) $@
