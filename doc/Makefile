LATEX = latex
DVIPDF = dvipdf
BIBTEX = bibtex
PDFLATEX = pdflatex

all: dating dp

dating: dating.tex
	$(LATEX) $+
	$(BIBTEX) $@
	$(DVIPDF) $@

dp: dp.tex
	$(LATEX) $+
	$(BIBTEX) $@
	$(DVIPDF) $@
