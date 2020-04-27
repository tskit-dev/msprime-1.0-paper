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

clean:
	make -C sims clean
	rm -fR data
