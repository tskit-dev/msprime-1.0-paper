DATA=data/arg.csv\
     data/mutations_perf.csv

PERF_FIGURES=figures/mutations-perf.pdf\
	figures/arg.pdf\
	figures/sweeps-perf.pdf\
	figures/gc-perf.pdf\
	figures/dtwf-perf.pdf

FIGURES=$(PERF_FIGURES)\
	figures/mutated_tree.pdf\
	figures/unmutated_tree.pdf\
	figures/example_tree_sequence.pdf

paper.pdf: paper.tex paper.bib ${DATA} ${FIGURES}
	pdflatex paper.tex
	bibtex paper
	pdflatex paper.tex
	pdflatex paper.tex

$(PERF_FIGURES):
	python3 evaluation/plot.py $*

figures/mutated_tree.svg figures/unmutated_tree.svg : pretty_pictures.py
	python3 $<

%.pdf : %.svg
	inkscape $< --export-pdf=$@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

data/arg.csv:
	python3 evaluation/generate_arg_data.py

data/mutations_perf.csv:
	python3 evaluation/generate_mutations_perf_data.py generate-trees
	python3 evaluation/generate_mutations_perf_data.py benchmark-on-trees

paper.ps: paper.dvi
	dvips paper

paper.dvi: paper.tex paper.bib
	latex paper.tex
	bibtex paper
	latex paper.tex
	latex paper.tex

clean:
	rm -f *.log *.dvi *.aux
	rm -f *.blg *.bbl
	rm -f *.eps *.[1-9]
	rm -f src/*.mpx *.mpx

mrproper: clean
	rm -f *.ps *.pdf
