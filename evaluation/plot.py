import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
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
    df = df.groupby(["sample_size", "tool"]).mean().reset_index()
    print(df)

    df.memory /= 1024 ** 3
    df.user_time /= 3600

    lines = {}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    for tool in ["SimBac", "fastSimBac", "msprime"]:
        dft = df[df.tool == tool]
        (line,) = ax1.plot(dft.sample_size, dft.user_time, label=tool)
        ax2.plot(dft.sample_size, dft.memory, label=tool, color=line.get_color())
        lines[tool] = line

    ax1.set_title("(A)")
    ax2.set_title("(B)")
    ax1.set_xlabel("Sample size")
    ax1.set_ylabel("Time (hours)")
    ax1.set_ylim(0, None)
    ax2.set_ylim(0, None)
    ax2.set_xlabel("Sample size")
    ax2.set_ylabel("Memory (GiB)")
    ax1.legend()

    dfmsp = df[df["tool"] == "msprime"]
    largest_n = np.array(dfmsp.sample_size)[-1]
    largest_value = np.array(dfmsp.user_time)[-1]
    ax1.plot([largest_n], [largest_value], "o", color=lines["msprime"].get_color())
    ax1.annotate(
        f"{round(largest_value * 60)} mins",
        textcoords="offset points",
        xytext=(-30, 15),
        xy=(largest_n, largest_value),
        xycoords="data",
    )
    plt.tight_layout()

    save("gc-perf")


@click.command()
def sweeps_perf():
    """
    Plot the sweeps benchmark.
    """
    df = pd.read_csv("data/sweeps_perf.csv")
    df = df.groupby(["L", "tool"]).mean().reset_index()
    print(df)

    df.memory /= 1024 ** 3
    df.L /= 1000
    # discoal has very high systime usage, so we need to
    # include it here as well.
    df["time"] = (df.user_time + df.sys_time) / 60

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    lines = {}
    for tool in ["msprime", "discoal"]:
        dft = df[df.tool == tool]
        (line,) = ax1.plot(dft.L, dft.time, label=tool)
        lines[tool] = line
        ax2.plot(dft.L, dft.memory, label=tool)

    ax1.set_title("(A)")
    ax2.set_title("(B)")
    ax1.set_xlabel("Sequence length (Kilobases)")
    ax1.set_ylabel("Time (minutes)")
    ax2.set_xlabel("Sequence length (Kilobases)")
    ax2.set_ylabel("Memory (GiB)")
    ax1.legend()

    dfmsp = df[df["tool"] == "msprime"]
    largest_L = np.array(dfmsp.L)[-1]
    largest_value = np.array(dfmsp.time)[-1]
    ax1.plot([largest_L], [largest_value], "o", color=lines["msprime"].get_color())
    ax1.annotate(
        f"{round(largest_value * 60)} seconds",
        textcoords="offset points",
        xytext=(-30, 15),
        xy=(largest_L, largest_value),
        xycoords="data",
    )
    largest_value = np.array(dfmsp.memory)[-1]
    ax2.plot([largest_L], [largest_value], "o", color=lines["msprime"].get_color())
    ax2.annotate(
        f"{round(largest_value * 1024)} MiB",
        textcoords="offset points",
        xytext=(-30, 15),
        xy=(largest_L, largest_value),
        xycoords="data",
    )
    plt.tight_layout()
    save("sweeps-perf")


@click.command()
def dtwf_perf():
    """
    Plot the DTWF benchark.
    """

    df = pd.read_csv("data/dtwf-perf.csv")
    df = df.groupby(["L", "tool"]).mean().reset_index()
    # print(df)

    df.memory /= 1024 ** 3
    label_map = {"ARGON": "ARGON", "msprime": "DTWF", "hybrid": "DTWF + Hudson"}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    for tool in ["ARGON", "msprime", "hybrid"]:
        dft = df[df.tool == tool]
        ax1.plot(dft.L, dft.user_time, label=label_map[tool])
        ax2.plot(dft.L, dft.memory, label=label_map[tool])

    ax1.set_xlabel("Sequence length (Megabases)")
    ax1.set_ylabel("Time (seconds)")
    ax2.set_xlabel("Sequence length (Megabases)")
    ax2.set_ylabel("Memory (GiB)")
    ax1.legend()
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
    print(f"Dromel rho {dromel_chr2l:.2g}")

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

            ax.set_xlabel("$N_e L$ (= scaled recombination rate $\\rho/4$)")

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
                # The quadratic is still fit to rho, so don't divide by 4 here
                fitted_quadratic(dromel_chr2l) / 3600,
                "hours",
            )

    axes[0].legend(handles=legend_adds)
    axes[1].legend()
    axes[0].set_title("A")
    axes[1].set_title("B")

    save("ancestry-perf")


cli.add_command(mutations_perf)
cli.add_command(gc_perf)
cli.add_command(sweeps_perf)
cli.add_command(dtwf_perf)
cli.add_command(ancestry_perf)

if __name__ == "__main__":
    cli()
