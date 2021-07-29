
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.optimize
import click


@click.group()
def cli():
    pass


def save(name):
    plt.tight_layout()
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
    # Get the mean among replicates
    df.time /= 3600
    dfg = df.groupby(["sample_size", "tool"]).mean().reset_index()

    fig, ax1 = plt.subplots(1, 1)
    lines = {}
    for tool in ["SimBac", "fastSimBac", "msprime"]:
        dft = dfg[dfg["tool"] == tool]
        (line,) = ax1.plot(dft.sample_size, dft.time, label=tool)
        lines[tool] = line

    ax1.set_xlabel("Sample size")
    ax1.set_ylabel("Time (hours)")
    # We could set a log-scale on the y-axis here but it really
    # doesn't make much difference

    dfmsp = dfg[dfg["tool"] == "msprime"]
    largest_n = np.array(dfmsp.sample_size)[-1]
    largest_value = np.array(dfmsp.time)[-1]
    ax1.plot([largest_n], [largest_value], "o", color=lines["msprime"].get_color())
    ax1.annotate(
        f"{round(largest_value * 60)} mins",
        textcoords="offset points",
        xytext=(-30, 15),
        xy=(largest_n, largest_value),
        xycoords="data",
    )
    ax1.legend()
    plt.tight_layout()
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
def ancestry_perf():
    """
    Plot the ancestry benchark.
    """

    df = pd.read_csv("data/ancestry-perf.csv", sep=",")

    def objective(X, a, b):
        r, n = X
        return a * (r ** 2) * (np.log(n)**2) + b

    # Constants aren't going to matter here since we're fitting, but not
    # to worry.
    rho = np.array(4 * df["N"] * df["L"] / 1e8)
    n = np.array(df["num_samples"])
    T = np.array(df["time"])
    popt, _ = scipy.optimize.curve_fit(objective, [rho, n], T, [1, 0])

    fig, axes = plt.subplots(1, 2, figsize=(6, 3))

    # Plot the times
    fig, axes = plt.subplots(1, 2, figsize=(6, 3))
    rgb = matplotlib.cm.get_cmap("Set1")(np.linspace(0.0, 1.0, len(set(df["N"]))))
    time_lims = [5, 1e12]
    for ax, tl  in zip(axes, time_lims):
        legend_adds = []
        for ns, m in zip((1000, 100000), ("o", "v")):
            legend_adds.append(
                    matplotlib.lines.Line2D(
                        [], [], color='grey', marker=m,
                        linestyle="none",
                        label=f"num samples={ns}"
                    )
            )
            ut = np.logical_and(
                    df["time"] <= tl,
                    df["num_samples"] == ns
            )
            # fit a qudratic with no intercept:
            # in R this would be lm(time ~ 0 + L : N + I((N*L)^2))
            rho = 4 * df["N"] * df["L"] / 1e8
            X = np.empty((sum(ut), 3))
            X[:,0] = rho[ut]
            X[:,1] = ( rho ** 2 )[ut]
            X[:,2] = 1
            b, _, _, _ = np.linalg.lstsq(X, df["time"][ut], rcond=None)

            ax.set_xlabel("rho (4NL/1e8)")
            ax.set_ylabel("time (seconds)")
            Nvals = sorted(list(set(df["N"])))
            for k, Nval in enumerate(Nvals):
                utN = np.logical_and(ut, df["N"] == Nval)
                if np.sum(utN) > 0:
                    sargs = {}
                    if m == "o":
                        sargs["label"] = f"N={Nval}"
                    ax.scatter(
                        rho[utN],
                        df["time"][utN],
                        color=rgb[k],
                        marker=m,
                        **sargs
                    )

            xx = np.linspace(0, 1.05 * max(X[:, 0]), 51)

            ax.plot(xx, [objective((x,ns), *popt) for x in xx], color="red")

            # Note the two quadratic curves are not the same!
            # pargs = {}
            # if m == "o":
            #     pargs["label"] = "quadratic"
            # ax.plot(xx, b[2] + b[0] * xx + b[1] * (xx ** 2), color="black", **pargs)
            print(f"Times less than {tl}: "
                  f"{b[2]:.2f} + {b[0]} * rho + {b[1]} * rho^2")

    axes[0].legend(handles=legend_adds, prop={'size': 6})
    axes[1].legend(prop={'size': 6})
    save("ancestry-perf")


@click.command()
def ancestry_perf_expected():
    """
    Compare the ancestry benchark to Hein, Schierup & Wiuf:
    recombination rate is 1e-8, so rho = N * L * 1e-8
    """

    df = pd.read_csv(
        "data/ancestry-perf.csv", sep=",",
    )

    h = np.cumsum(1/np.arange(1, np.max(df["num_samples"]) + 1))
    rho = df["L"] * df["N"] * 1e-8
    expected = rho * (rho + 1) * (h[df["num_samples"].astype("int") - 1] ** 2)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(6, 3))

    ax0.set_xlabel("expected time (theory)")
    ax0.set_ylabel("observed time (seconds)")
    for ns, m in zip((1000, 100000), ("o", "v")):
        ut = (df["num_samples"] == ns)
        ax0.scatter(
                expected[ut], df["time"][ut],
                marker=m, label=f"num_samples={ns}",
        )
    ax0.legend(prop={'size': 6})

    ax1.set_xlabel("rho (4NL/1e8)")
    ax1.set_ylabel("observed / expected time")
    for ns, m in zip((1000, 100000), ("o", "v")):
        ut = (df["num_samples"] == ns)
        ax1.scatter(
                (4 * df["N"] * df["L"] / 1e8)[ut],
                (df["time"] / expected)[ut],
                marker=m,
        )
    ax1.set_yscale("log")

    save("ancestry-perf-expected")


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
cli.add_command(ancestry_perf)
cli.add_command(ancestry_perf_expected)
cli.add_command(arg)

if __name__ == "__main__":
    cli()
