This directory contains scripts to automatically install
and benchmark existing simulation tools against msprime. These
tools include:

- ms
- msms
- scrm
- discoal
- ARGON

To install and build the simulators, run:
```
make sims
```
This also installs the version of msprime specified in
`requirements.txt`

The `generate_performance_data.py` script runs comparisons of the
standard coalescent, DTWF, SMC approximation, and selective sweeps.
To generate performance data for a particular model, run
```
python generate_performance_data.py <model>
```
To generate data for all models, run
```
make data
```

To plot the data from all simulations, run
```
make plots
```
