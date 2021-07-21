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


def save(name):
    plt.savefig(f"figures/{name}.png")
    plt.savefig(f"figures/{name}.pdf")


@click.command()
def mutations_perf():
    """
    Plot the mutations benchmark.
    """
    df = pd.read_csv("data/mutations_perf.csv")
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
    save("mutations-perf")


@click.command()
def gc_perf():
    """
    Plot the gene-conversion benchmark.
    """
    df = pd.read_csv("data/gc-perf.csv")

    fig, ax1 = plt.subplots(1, 1)
    for tool in set(df.tool):
        dft = df[df.tool == tool]
        ax1.plot(dft.sample_size, dft.time, label=tool)

    ax1.set_xlabel("Sample size")
    ax1.set_ylabel("Time (seconds)")

    ax1.legend()
    save("gc-perf")


@click.command()
def sweeps_perf():
    """
    Plot the sweeps benchmark.
    """
    df = pd.read_csv("data/sweeps_perf.csv")

    fig, ax1 = plt.subplots(1, 1)
    for tool in set(df.tool):
        dft = df[df.tool == tool]
        ax1.plot(dft.L, dft.time, label=tool)

    ax1.set_xlabel("Sequence length")
    ax1.set_ylabel("Time (seconds)")

    ax1.legend()
    save("sweeps-perf")


@click.command()
def dtwf_perf():
    """
    Plot the DTWF benchark.
    """

    df = pd.read_csv(
        "data/dtwf-perf.csv", sep="\t", usecols=["length", "program", "runtime"]
    )

    # Make a new dataframe with just the runtime means.
    means = df.groupby(["length", "program"]).mean().reset_index()
    sample_sizes = [0, 1, 2, 3, 4]
    sample_size_labs = [
        r"$10^{-2}$",
        r"$10^{-1}$",
        r"$10^{0}$",
        r"$10^{1}$",
        r"$10^{2}$",
    ]

    # print(means)

    # # Plot the means
    plt.plot(
        sample_sizes,
        np.log10(means[means["program"] == "msprime-DTWF"]["runtime"]),
        label="msprime-DTWF",
        marker="o",
        color="C1",
    )
    plt.plot(
        sample_sizes,
        np.log10(means[means["program"] == "msprime-hybrid"]["runtime"]),
        label="msprime-hybrid",
        marker="o",
        color="C1",
        linestyle="dashed",
    )
    plt.plot(
        sample_sizes,
        np.log10(means[means["program"] == "argon"]["runtime"]),
        label="argon",
        marker="o",
        color="C0",
    )
    plt.legend(loc="upper left")
    plt.xlabel("Chromosome lengths (Mb)")
    plt.ylabel("Runtimes (secs)")
    plt.title("Simulation runtimes for different genome lengths")
    plt.xticks(sample_sizes, labels=sample_size_labs)
    plt.yticks(
        [-2, -1, 0, 1, 2],
        labels=[r"$10^{-2}$", r"$10^{1}$", r"$10^{0}$", r"$10^{1}$", r"$10^{2}$"],
    )
    plt.grid(b=True, linewidth=0.5, linestyle="dotted")
    save("dtwf-perf")


@click.command()
def arg():
    """
    Plot the ARG size illustration figure.
    """
    df = pd.read_csv("data/arg.csv")

    fig, ax1 = plt.subplots(1, 1)
    ax1.plot(df.L, df.arg_nodes, label="Fraction of ARG nodes")
    ax1.plot(df.L, df.size_ratio, label="size(tree sequence) / size(ARG)")
    ax1.set_xlabel("Sequence length")
    ax1.legend()
    save("arg")


cli.add_command(mutations_perf)
cli.add_command(gc_perf)
cli.add_command(sweeps_perf)
cli.add_command(dtwf_perf)
cli.add_command(arg)

if __name__ == "__main__":
    cli()
