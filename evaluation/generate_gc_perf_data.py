import concurrent.futures
import subprocess
import tempfile
import time

import cpuinfo
import pandas as pd
import numpy as np
import click
import msprime
import matplotlib.pyplot as plt
import statsmodels.api as sm


def run_simbac(*, sample_size, L, gc_rate, gc_tract_length, count_trees=False):

    # using R=2*gc_rate as gene conversion/recombination rate as SimBac uses R/2
    R = gc_rate * 2
    # Set theta to 0 to remove mutations (defaults to 0.01)
    args = f"-N {sample_size} -B {int(L)} -R {R} -D {gc_tract_length} -T 0"
    if count_trees:
        # This is cleaned up automatically when the variable goes out of scope
        tempdir = tempfile.TemporaryDirectory()
        outfile = tempdir.name + "/out.newick"
        args += " -l " + outfile

    # print("running", args)
    subprocess.run(
        "./tools/simbac " + args, shell=True, stdout=subprocess.DEVNULL, check=True
    )
    num_trees = 0
    if count_trees:
        with open(outfile) as f:
            for line in f:
                num_trees += 1

    return num_trees


def run_fastsimbac(
    *, sample_size, L, gc_rate, gc_tract_length, set_seed=0, count_trees=False
):

    # using R=2*gc_rate as gene conversion/recombination rate
    R = gc_rate * 2
    # Set theta to 0 to remove mutations (defaults to 0.01)
    args = f"{sample_size} {int(L)} -r {R} {gc_tract_length} -t 0"
    if count_trees:
        args += " -T -b 500"
    if set_seed > 0:
        args += " -s " + str(set_seed)

    # print("running", args)
    output = subprocess.run(
        "./tools/fastSimBac " + args, shell=True, check=True, capture_output=True
    )
    # print(output.stdout)
    num_trees = 0
    if count_trees:
        for line in output.stdout.splitlines():
            if line.startswith(b"["):
                num_trees += 1

    return num_trees


def run_msprime(*, sample_size, L, gc_rate, gc_tract_length, ret_breakpoints=True):
    sim = msprime.sim_ancestry(
        samples=sample_size,
        sequence_length=L,
        ploidy=1,
        gene_conversion_rate=gc_rate,
        gene_conversion_tract_length=gc_tract_length,
    )
    treenumber = sim.num_trees

    # We use an internal msprime API here because we want to get at the
    # number of breakpoints, not the distinct trees.
    if ret_breakpoints:
        sim = msprime.ancestry._parse_sim_ancestry(
            samples=sample_size,
            sequence_length=L,
            ploidy=1,
            gene_conversion_rate=gc_rate,
            gene_conversion_tract_length=gc_tract_length,
        )
        sim.run()
        breakpointnumber = sim.num_breakpoints
        return treenumber, breakpointnumber

    return treenumber


def run_benchmark(work):
    tool, sample_size = work

    L = 4_500_000
    gc_rate = 0.015
    gc_tract_length = 500

    func = tool_map[tool]
    before = time.perf_counter()
    func(
        sample_size=sample_size, L=L, gc_rate=gc_rate, gc_tract_length=gc_tract_length,
    )
    duration = time.perf_counter() - before
    return {
        "sample_size": sample_size,
        "tool": tool,
        "time": duration,
    }


@click.command()
@click.option("--replicates", type=int, default=1000)
@click.option("--sample-size", type=int, default=10)
def validate(replicates, sample_size):
    """
    Validate that we are simulating the same things in the simulators
    by running some replicates and plotting the distributions of the
    number of output trees.
    """
    L = 1000
    gc_rate = 0.015
    gc_tract_length = 10

    nt_simbac = np.zeros(replicates)
    nt_fastsimbac = np.zeros(replicates)
    nt_msprime = np.zeros(replicates)
    nb_msprime = np.zeros(replicates)

    with click.progressbar(range(replicates)) as bar:
        for j in bar:
            nt_simbac[j] = run_simbac(
                sample_size=sample_size,
                L=L,
                gc_rate=gc_rate,
                gc_tract_length=gc_tract_length,
                count_trees=True,
            )
            nt_fastsimbac[j] = run_fastsimbac(
                sample_size=sample_size,
                L=L,
                gc_rate=gc_rate,
                gc_tract_length=gc_tract_length,
                set_seed=j,
                count_trees=True,
            )
            nt_msprime[j], nb_msprime[j] = run_msprime(
                sample_size=sample_size,
                L=L,
                gc_rate=gc_rate,
                gc_tract_length=gc_tract_length,
                ret_breakpoints=True,
            )
    print(
        "mean number of trees:",
        "simbac=",
        np.mean(nt_simbac),
        "fastsimbac=",
        np.mean(nt_fastsimbac),
        "msprime trees=",
        np.mean(nt_msprime),
        "msprime breakpoints=",
        np.mean(nb_msprime),
    )

    sm.graphics.qqplot(nt_simbac)
    sm.qqplot_2samples(nt_simbac, nb_msprime, line="45")
    plt.xlabel("simbac")
    plt.ylabel("msprime")
    plt.savefig("figures/verify_simbac_v_msprime.png")

    plt.close("all")

    sm.graphics.qqplot(nt_fastsimbac)
    sm.qqplot_2samples(nt_fastsimbac, nt_msprime, line="45")
    plt.xlabel("fastsimbac")
    plt.ylabel("msprime")
    plt.savefig("figures/verify_fastsimbac_v_msprime.png")


tool_map = {
    "msprime": run_msprime,
    "SimBac": run_simbac,
    "fastSimBac": run_fastsimbac,
}


@click.command()
@click.option("--replicates", type=int, default=5)
@click.option("--processes", type=int, default=None)
def benchmark_ecoli(replicates, processes):
    """
    Runs the benchmarks for E-coli simulations.
    """

    cpu = cpuinfo.get_cpu_info()
    with open("data/gc_perf_cpu.txt", "w") as f:
        for k, v in cpu.items():
            print(k, "\t", v, file=f)

    sample_sizes = np.linspace(10, 500, 20).astype(int)
    work = []
    for sample_size in sample_sizes:
        for name in tool_map.keys():
            if sample_size > 50 and name == "SimBac":
                continue
            work.extend([(name, sample_size)] * replicates)
    data = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=processes) as executor:
        futures = [executor.submit(run_benchmark, item) for item in work]
        for future in concurrent.futures.as_completed(futures):
            data.append(future.result())
            print(data[-1])
            df = pd.DataFrame(data)
            df.to_csv("data/gc-perf.csv")


@click.group()
def cli():
    pass


cli.add_command(validate)
cli.add_command(benchmark_ecoli)

if __name__ == "__main__":
    cli()
