cover_letter3.pdf: cover_letter3.tex
	pdflatex cover_letter3.tex

cover_letter2.pdf: cover_letter2.tex
	pdflatex cover_letter2.tex

cover_letter.pdf: cover_letter.tex
	pdflatex cover_letter.tex

clean:
	rm -f *.log *.dvi *.aux
	rm -f *.blg *.bbl
	rm -f src/*.mpx *.mpx

mrproper: clean
	rm -f *.ps *.pdf
