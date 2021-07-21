#!/usr/bin/env python3

import tskit, msprime
import subprocess
import pandas as pd
import time

# Methods for simulating with ARGON and msprime.

def sim_argon(sample_size, chrom_length, min_length=0, no_rec=False,
    classpath="/Users/gtsambos/Software/ARGON/", ibd=True, trees=True):
    # Defaults: mutation rate 1.65e-8, rec rate 1e-8
    diploid_size = 2*sample_size
    arguments = ['java',
            '-jar', classpath + 'ARGON.0.1.jar',
            '-N', '10000',
            '-pop', '1', '%i' % diploid_size,
            '-size', '%i' % int(chrom_length/1000000),
            '-mut', '0',
            '-quiet'
            ]
    if trees:
        arguments.extend(['-trees'])
    if ibd:
        arguments.extend(['-IBD', '%i' % min_length])
    if no_rec:
        arguments.extend(['-rec', '0'])

    process = subprocess.run(arguments)
    stdout, stderr = process.communicate()
    return (stdout, stderr)

# Methods for simulating with msprime.

def sim_msprime(sample_size, chrom_length, min_length=0, max_time=None, ibd=True,
    no_rec=False):
    if no_rec:
        rho = 0
    else:
        rho = 1e-8
    ts = msprime.sim_ancestry(
        samples=sample_size, sequence_length=chrom_length,
        discrete_genome=True,
        recombination_rate=rho, population_size=5000,
        model='dtwf')
    # ts = msprime.sim_mutations(ts, rate=1e-20)
    if not ibd:
        return ts
    else:
        samples = list(itertools.combinations(range(sample_size), 2))
        ibd = ts.tables.find_ibd(
            min_length=min_length, samples=samples, max_time=max_time)
        # print("Number of simulated trees:", ts.num_trees)
        return ibd

def sim_msprime_hybrid(sample_size, chrom_length, min_length=0, max_time=None, ibd=True,
    no_rec=False):
    if no_rec:
        rho = 0
    else:
        rho = 1e-8
    ts = msprime.sim_ancestry(
        samples=sample_size, sequence_length=chrom_length,
        discrete_genome=True,
        recombination_rate=rho, population_size=5000,
        model = [msprime.DiscreteTimeWrightFisher(duration=100),
        msprime.StandardCoalescent()])
    # ts = msprime.sim_mutations(ts, rate=1e-20)
    if not ibd:
        return ts
    else:
        samples = list(itertools.combinations(range(sample_size), 2))
        ibd = ts.tables.find_ibd(
            min_length=min_length, samples=samples, max_time=max_time)
        # print("Number of simulated trees:", ts.num_trees)
        return ibd

# Define constants
num_replicates = 20
sample_size = 400
# chrom_lengths = [500000, 5000000, 50000000, 500000000]
chrom_lengths = [10000, 100000, 1000000, 10000000, 100000000]
cutoff=0

# Initialise lists to turn into outputted dataframe.
df_samples = []
df_length = []
df_cutoff = []
df_program = []
df_runtime = []

# Loop over sample sizes
for c in chrom_lengths:

    # Loop over replicates
    for r in range(num_replicates):

        ### msprime DTWF

        time_start = time.time()
        utils.sim_msprime(sample_size, c, min_length=0, ibd=False)
        time_end = time.time()

        # Append to lists.
        df_length.append(c)
        df_program.append('msprime-DTWF')
        df_runtime.append(time_end-time_start)

        ### msprime hybrid: 100 gens DTWF, rest Hudson

        time_start = time.time()
        utils.sim_msprime_hybrid(sample_size, c, min_length=0, ibd=False)
        time_end = time.time()

        # Append to lists.
        df_length.append(c)
        df_program.append('msprime-hybrid')
        df_runtime.append(time_end-time_start)

        #### ARGON

        time_start = time.time()
        stdo, stde = utils.sim_argon(sample_size, c, min_length=0, ibd=False, trees=False)
        time_end = time.time()

        # Append to lists.
        df_length.append(c)
        df_program.append('argon')
        df_runtime.append(time_end-time_start)

# Make pandas dataframe and save it.
out = pd.DataFrame({
    'length' : df_length,
    'program' : df_program,
    'runtime' : df_runtime
    })
print(out)

out.to_csv('output/x.msprime-paper-simulate-genomelengths.csv', sep="\t")

