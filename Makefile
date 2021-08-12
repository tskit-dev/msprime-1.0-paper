DATA=data/mutations_perf.csv

FIGURES=figures/mutations-perf.pdf\
	figures/sweeps-perf.pdf\
	figures/gc-perf.pdf\
	figures/dtwf-perf.pdf\
	figures/ancestry-perf.pdf

ILLUSTRATIONS=\
	illustrations/mutated_tree.pdf\
	illustrations/example_tree_sequence.pdf

paper.pdf: paper.tex paper.bib ${DATA} ${FIGURES} ${ILLUSTRATIONS}
	pdflatex paper.tex
	bibtex paper
	pdflatex paper.tex
	pdflatex paper.tex

figures/%.pdf:
	python3 evaluation/plot.py $*

illustrations/mutated_tree.svg: pretty_pictures.py
	python3 $< mutated-tree

illustrations/arg-ts.svg: pretty_pictures.py
	python3 $< arg-ts

%.pdf : %.svg
	# Have to use Chromium here to get rid of PDF pages
	chromium --headless --no-margins --print-to-pdf=$@ $<
	# inkscape $< --export-pdf=$@

%.pdf : %.ink.svg
	inkscape $< --export-pdf=$@

data/mutations_perf.csv:
	python3 evaluation/generate_mutations_perf_data.py generate-trees
	python3 evaluation/generate_mutations_perf_data.py benchmark-on-trees

data/ancestry_perf.csv:
	python3 evaluation/generation_ancestry_perf_data.py

paper.ps: paper.dvi
	dvips paper

paper.dvi: paper.tex paper.bib
	latex paper.tex
	bibtex paper
	latex paper.tex
	latex paper.tex

.PHONY: spellcheck
spellcheck: aspell.conf
	aspell --conf ./aspell.conf --check paper.tex

clean:
	rm -f *.log *.dvi *.aux
	rm -f *.blg *.bbl
	rm -f *.eps *.[1-9]
	rm -f src/*.mpx *.mpx

mrproper: clean
	rm -f *.ps *.pdf
