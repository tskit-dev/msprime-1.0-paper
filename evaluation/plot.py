import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import click

# Main text size is 9pt
plt.rcParams.update({"font.size": 7})
plt.rcParams.update({"legend.fontsize": 6})
plt.rcParams.update({"lines.markersize": 4})


@click.group()
def cli():
    pass


def save(name):
    plt.tight_layout()
    plt.savefig(f"figures/{name}.png")
    plt.savefig(f"figures/{name}.pdf")


def two_panel_fig(**kwargs):
    # The columnwidth of the genetics format is ~250pt, which is
    # 3 15/32 inch, = 3.46
    width = 3.46
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(width, width / 2), **kwargs)
    ax1.set_title("(A)")
    ax2.set_title("(B)")
    return fig, (ax1, ax2)

def two_panel_fig_single_col(**kwargs):
    width = 6
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(width, width / 2), **kwargs)
    ax1.set_title("(A)")
    ax2.set_title("(B)")
    return fig, (ax1, ax2)


@click.command()
def mutations_perf():
    """
    Plot the mutations benchmark.
    """
    df = pd.read_csv("data/mutations_perf.csv")
    df.L /= 1e6
    rates = np.unique(df.rate)

    fig, (ax1, ax2) = two_panel_fig(sharey=True)
    ax1.set_xlabel("Sample size (haploid)")
    dfL = df[df.L == 10]
    markers= ["+", "x", "1"]
    for marker, rate in zip(markers, rates):
        dfr = dfL[dfL.rate == rate].sort_values("n")
        rate_exp = int(np.log10(rate))
        ax1.plot(dfr.n, dfr.time, marker=marker,
        label=f"Mutation rate=$10^{{{rate_exp}}}$")

    ax2.set_xlabel("Sequence length (Megabases)")
    ax1.set_ylabel("Time (seconds)")

    dfn = df[df.n == 1000].sort_values("L")
    for marker, rate in zip(markers, rates):
        dfr = dfn[dfn.rate == rate]
        ax2.plot(dfr.L, dfr.time, marker=marker, label=f"Mutation rate={rate}")
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
    fig, (ax1, ax2) = two_panel_fig()
    marker_map = {
        "msprime": ".",
        "SimBac": "+",
        "fastSimBac": "x",
    }

    for tool in ["msprime", "SimBac", "fastSimBac"]:
        dft = df[df.tool == tool]
        (line,) = ax1.plot(dft.sample_size, dft.user_time,
                marker=marker_map[tool], label=tool)
        ax2.plot(dft.sample_size, dft.memory, label=tool,
                marker=marker_map[tool], color=line.get_color())
        lines[tool] = line

    ax1.set_xlabel("Sample size (haploid)")
    ax1.set_ylabel("Time (hours)")
    ax2.set_xlabel("Sample size (haploid)")
    ax2.set_ylabel("Memory (GiB)")

    dfmsp = df[df["tool"] == "msprime"]
    largest_n = np.array(dfmsp.sample_size)[-1]
    largest_value = np.array(dfmsp.user_time)[-1]
    ax1.plot([largest_n], [largest_value], marker="*", color="black")
    ax1.annotate(
        f"{round(largest_value * 60)} mins",
        textcoords="offset points",
        xytext=(-25, 5),
        xy=(largest_n, largest_value),
        xycoords="data",
    )
    ax2.set_ylim(bottom=0)
    ax2.legend()

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
    df.L /= 1000**2
    # discoal has very high systime usage, so we need to
    # include it here as well.
    df["time"] = (df.user_time + df.sys_time) / 60

    fig, (ax1, ax2) = two_panel_fig()
    lines = {}
    marker_map = {
        "msprime": ".",
        "discoal": "x",
    }

    for tool in ["msprime", "discoal"]:
        dft = df[df.tool == tool]
        (line,) = ax1.plot(dft.L, dft.time, marker=marker_map[tool], label=tool)
        lines[tool] = line
        ax2.plot(dft.L, dft.memory, marker=marker_map[tool], label=tool)

    ax1.set_xlabel("Sequence length (Megabases)")
    ax1.set_ylabel("Time (minutes)")
    ax2.set_xlabel("Sequence length (Megabases)")
    ax2.set_ylabel("Memory (GiB)")
    ax1.legend()

    dfmsp = df[df["tool"] == "msprime"]
    largest_L = np.array(dfmsp.L)[-1]
    largest_value = np.array(dfmsp.time)[-1]
    ax1.plot([largest_L], [largest_value], "*", color="black")
    ax1.annotate(
        f"{largest_value * 60:.1f} seconds",
        textcoords="offset points",
        xytext=(-39, 5),
        xy=(largest_L, largest_value),
        xycoords="data",
    )
    largest_value = np.array(dfmsp.memory)[-1]
    ax2.plot([largest_L], [largest_value], "*", color="black")
    ax2.annotate(
        f"{round(largest_value * 1024)} MiB",
        textcoords="offset points",
        xytext=(-21, 5),
        xy=(largest_L, largest_value),
        xycoords="data",
    )
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
    label_map = {"ARGON": "ARGON", "msprime": "DTWF", "hybrid": "Hybrid"}

    fig, (ax1, ax2) = two_panel_fig()
    lines = {}
    marker_map = {"msprime": ".",
            "hybrid": "+", "ARGON": "x"}

    for tool in ["msprime", "hybrid", "ARGON"]:
        dft = df[df.tool == tool]
        ax1.plot(dft.L, dft.user_time, marker=marker_map[tool], label=label_map[tool])
        (line,) = ax2.plot(dft.L, dft.memory, marker=marker_map[tool], label=label_map[tool])
        lines[tool] = line

    ax1.set_xlabel("Sequence length (Megabases)")
    ax1.set_ylabel("Time (seconds)")
    ax2.set_xlabel("Sequence length (Megabases)")
    ax2.set_ylabel("Memory (GiB)")

    dfmsp = df[df["tool"] == "hybrid"]
    largest_L = np.array(dfmsp.L)[-1]
    largest_value = np.array(dfmsp.memory)[-1]
    ax2.plot([largest_L], [largest_value], "*", color="black")
    ax2.annotate(
        f"{round(largest_value * 1024)} MiB",
        textcoords="offset points",
        xytext=(-26, 5),
        xy=(largest_L, largest_value),
        xycoords="data",
    )
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

    fig, axes = two_panel_fig_single_col()

    def annotate_rho(ax, rho, x_offset, y_offset, text):
        ax.axvline(rho / 4, color="0.8", linestyle="-.")
        ax.text(rho / 4 + x_offset, y_offset, text, fontstyle="italic")

    annotate_rho(axes[0], aratha_chr1, 100, 1.5, "Arabidopsis\nthaliana")
    annotate_rho(axes[1], human_chr1, -5000, 100, "Homo sapiens")
    annotate_rho(axes[1], canfam_chr1, 1000, 0, "Canis familiaris")
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

            # xx = np.linspace(0, 1.05 * max(X[:, 0]), 51)
            xx = np.linspace(0, max(X[:, 0]), 51)

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
        max_x = np.max(rho[df["time"] <= tl] / 4)
        ticks = np.arange(
                0,
                1.1 * max_x,
                10 ** np.floor(np.log10(max_x)),
        )
        ax.set_xticks(ticks)
        ax.set_xticklabels(
                [f"{x:.0f}"for x in ticks]
        )

    prop = {"size": 7}
    axes[0].legend(handles=legend_adds, prop=prop)
    axes[1].legend(prop=prop)

    save("ancestry-perf")


cli.add_command(mutations_perf)
cli.add_command(gc_perf)
cli.add_command(sweeps_perf)
cli.add_command(dtwf_perf)

with matplotlib.rc_context({"font.size": 7}):
    cli.add_command(ancestry_perf)

if __name__ == "__main__":
    cli()
