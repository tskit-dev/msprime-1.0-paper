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

    # GRCh38, via stdpopsim
    human_rho = lambda L: int(4 * 10 ** 4 * L * 1e-8)
    human_chr1 = human_rho(248956422)

    # TAIR10, via stdpopsim
    r = 8.06e-10
    aratha_rho = lambda L: int(4 * 10 ** 4 * L * r)
    aratha_chr1 = aratha_rho(30427671)

    canfam_chr1 = 4 * 13000 * 122678785 * 7.636001498077e-09
    dromel_chr2l = 4 * 1720600 * 23513712 * 2.40462600791e-08

    # This is the figsize used in other two-panel plots
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    def annotate_rho(ax, rho, x_offset, text):
        ax.axvline(rho / 4, color="0.8", linestyle="-.")
        ax.text(rho / 4 + x_offset, 0, text, fontstyle="italic")

    annotate_rho(axes[0], aratha_chr1, 200, "Arabidopsis thaliana")
    annotate_rho(axes[1], human_chr1, 1000, "Homo sapiens")
    annotate_rho(axes[1], canfam_chr1, 1000, "Canis familiaris")
    # axes[1].axvline(dromel_chr2l)

    rgb = matplotlib.cm.get_cmap("Set1")(np.linspace(0.0, 1.0, len(set(df["N"]))))
    time_lims = [5, 1e12]
    for ax, tl in zip(axes, time_lims):
        legend_adds = [
            matplotlib.lines.Line2D(
                [],
                [],
                color="black",
                linestyle="-",
                label=f"quadratic",
            )
        ]
        for ns, m in zip((1000, 100000), ("o", "v")):
            legend_adds.append(
                matplotlib.lines.Line2D(
                    [],
                    [],
                    color="grey",
                    marker=m,
                    linestyle="none",
                    label=f"n={ns}",
                )
            )
            ut = np.logical_and(df["time"] <= tl, df["num_samples"] == ns)
            # fit a qudratic with no intercept:
            # in R this would be lm(time ~ 0 + L : N + I((N*L)^2))
            rho = 4 * df["N"] * df["L"] / 1e8
            X = np.empty((sum(ut), 3))
            X[:, 0] = rho[ut]
            X[:, 1] = (rho ** 2)[ut]
            X[:, 2] = 1
            b, _, _, _ = np.linalg.lstsq(X, df["time"][ut], rcond=None)

            ax.set_xlabel("Scaled recombination rate $\\rho/4 = N_e L$")

            def fitted_quadratic(x):
                return b[2] + b[0] * x + b[1] * (x ** 2)

            ax.set_ylabel("Time (seconds)")
            Nvals = sorted(list(set(df["N"])))
            for k, Nval in enumerate(Nvals):
                utN = np.logical_and(ut, df["N"] == Nval)
                if np.sum(utN) > 0:
                    sargs = {}
                    if m == "o":
                        sargs["label"] = f"$N_e={Nval}$"
                    ax.scatter(
                        rho[utN] / 4, df["time"][utN], color=rgb[k], marker=m, **sargs
                    )

            xx = np.linspace(0, 1.05 * max(X[:, 0]), 51)

            # Note the two quadratic curves are not the same!
            pargs = {}
            if m == "o":
                pargs["label"] = "quadratic"
            ax.plot(xx / 4, fitted_quadratic(xx), color="black", **pargs)
            # print(
            #     f"Times less than {tl}: " f"{b[2]:.2f} + {b[0]} * rho + {b[1]} * rho^2"
            # )
            print(
                "Predicted time for DroMel chr2L with n =",
                ns,
                "=",
                fitted_quadratic(dromel_chr2l / 4) / 3600,
                "hours",
            )

    axes[0].legend(handles=legend_adds)
    axes[1].legend()
    axes[0].set_title("A")
    axes[1].set_title("B")

    save("ancestry-perf")


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
cli.add_command(arg)

if __name__ == "__main__":
    cli()
