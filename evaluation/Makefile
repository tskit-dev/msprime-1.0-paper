all: sims data plots

.PHONY: sims data plots

sims:
	pip install -r requirements.txt
	make -C sims

data:
	mkdir -p data
	python3 generate_performance_data.py

plots:
	python3 plot_data.py data

# TODO clean these targets up depending on what we end up doing
# (or not) as the comparison benchmarks for mutation generation.

data/L_10kb_samples_10.trees:
	msp ancestry 10 -L 10000 --ploidy=1 -s 42 -o $@

data/L_10kb_samples_10k.trees:
	msp ancestry 10000 -L 10000 --ploidy=1 -s 42 -o $@

data/L_10kb_samples_100k.trees:
	msp ancestry 100000 -L 10000 --ploidy=1 -s 42 -o $@

data/L_10kb_samples_1m.trees:
	msp ancestry 1000000 -L 10000 --ploidy=1 -s 42 -o $@

data/L_30kb_samples_1m.trees:
	msp ancestry 1000000 -L 30000 --ploidy=1 -s 42 -o $@

# Benchmark seq-gen:
# Don't know how long it would take.
# ./sims/seq-gen -l10000 -m HKY data/L_1kb_samples_1m.nwk > /dev/null

clean:
	make -C sims clean
	rm -fR data
