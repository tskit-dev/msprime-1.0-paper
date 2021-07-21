#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['text.usetex'] = True

df = pd.read_csv("output/x.msprime-paper-simulate-genomelengths.csv", sep="\t",
    usecols=['length', 'program', 'runtime'])

# print(df)

# Make a new dataframe with just the runtime means.
means = df.groupby(['length', 'program']).mean().reset_index()
sample_sizes = [0, 1, 2, 3, 4]
sample_size_labs = [r'$10^{-2}$', r'$10^{-1}$', r'$10^{0}$', r'$10^{1}$', r'$10^{2}$']

# print(means)

# # Plot the means
plt.plot(sample_sizes, np.log10(means[means['program']=='msprime-DTWF']['runtime']),
    label='msprime-DTWF', marker='o', color='C1')
plt.plot(sample_sizes, np.log10(means[means['program']=='msprime-hybrid']['runtime']),
    label='msprime-hybrid', marker='o', color='C1', linestyle='dashed')
plt.plot(sample_sizes, np.log10(means[means['program']=='argon']['runtime']),
    label='argon', marker='o', color='C0')
plt.legend(loc='upper left')
plt.xlabel("Chromosome lengths (Mb)")
plt.ylabel("Runtimes (secs)")
plt.title("Simulation runtimes for different genome lengths")
plt.xticks(sample_sizes, labels=sample_size_labs)
plt.yticks([-2, -1, 0, 1, 2],
    labels=[r'$10^{-2}$', r'$10^{1}$', r'$10^{0}$', r'$10^{1}$', r'$10^{2}$'])
plt.grid(b=True, linewidth=0.5, linestyle='dotted')
plt.show()