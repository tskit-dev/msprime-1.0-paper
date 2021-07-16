import time
import pathlib

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import click
import pyvolve
import tskit
import msprime
import cpuinfo


@click.group()
def cli():
    pass


@click.command()
def mutation_perf():
    """
    Plot the mutations benchmark.
    """
    df = pd.read_csv("data/mutations.csv")
    df.L /= 1e6
    rates = np.unique(df.rate)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
    ax1.set_title("(A)")
    ax1.set_xlabel("Sample size")
    dfL = df[df.L == 10]
    for rate in rates:
        dfr = dfL[dfL.rate == rate].sort_values("n")
        ax1.plot(dfr.n, dfr.time, label=f"rate={rate}")

    ax2.set_title("(B)")
    ax2.set_xlabel("Sequence length (Megabases)")
    ax1.set_ylabel("Time (seconds)")

    dfn = df[df.n == 1000].sort_values("L")
    for rate in rates:
        dfr = dfn[dfn.rate == rate]
        ax2.plot(dfr.L, dfr.time, label=f"rate={rate}")
    print(dfn[dfn.L == 100])

    ax1.legend()
    plt.savefig("figures/mutation-perf.png")
    plt.savefig("figures/mutation-perf.pdf")


@click.command()
def gc_perf():
    """
    Plot the gene-conversion benchmark.
    """
    df = pd.read_csv("data/gc.csv")

    fig, ax1 = plt.subplots(1, 1)
    for tool in set(df.tool):
        dft = df[df.tool == tool]
        ax1.plot(dft.sample_size, dft.time, label=tool)

    ax1.set_xlabel("Sample size")
    ax1.set_ylabel("Time (seconds)")

    ax1.legend()
    plt.savefig("figures/gc-perf.png")
    plt.savefig("figures/gc-perf.pdf")


cli.add_command(mutation_perf)
cli.add_command(gc_perf)

if __name__ == "__main__":
    cli()
