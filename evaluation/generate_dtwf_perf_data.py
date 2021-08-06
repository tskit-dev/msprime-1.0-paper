import time
import subprocess
import concurrent.futures
import multiprocessing
import tempfile
import os
import resource

import matplotlib.pyplot as plt
import statsmodels.api as sm
import numpy as np
import click
import tskit
import msprime
import pandas as pd
import cpuinfo
import sys


def sim_argon(sample_size, chrom_length_mb, count_trees=False):
    tmpdir = tempfile.TemporaryDirectory()
    prefix = tmpdir.name + "/x"
    # Defaults: mutation rate 1.65e-8, rec rate 1e-8
    diploid_size = 2 * sample_size
    arguments = [
        "java",
        "-jar",
        "tools/ARGON.jar",
        "-N",
        "10000",
        "-pop",
        "1",
        f"{diploid_size}",
        "-size",
        f"{chrom_length_mb}",
        "-mut",
        "0",
        "-quiet",
        "-out",
        prefix,
    ]
    if count_trees:
        arguments.extend(["-trees"])
    # print("  ".join(arguments))
    subprocess.run(arguments, check=True, capture_output=True)
    num_trees = 0
    if count_trees:
        treefile = prefix + ".trees"
        # count the unique trees in the output. Technically we only want to
        # count *adjacent* identical trees, though.
        cmd = f"cut -f 4 {treefile} | uniq | wc -l"
        output = subprocess.run(cmd, shell=True, check=True, capture_output=True)
        num_trees = int(output.stdout)
    return num_trees


def sim_msprime(sample_size, chrom_length_mb, model="dtwf"):
    ts = msprime.sim_ancestry(
        samples=sample_size,
        sequence_length=chrom_length_mb * 10 ** 6,
        discrete_genome=True,
        recombination_rate=1e-8,
        population_size=5000,
        model=model,
    )
    return ts.num_trees


def sim_msprime_hybrid(sample_size, chrom_length_mb):
    return sim_msprime(
        sample_size,
        chrom_length_mb,
        model=[
            msprime.DiscreteTimeWrightFisher(duration=100),
            msprime.StandardCoalescent(),
        ],
    )


def process_resources():
    utime = 0
    stime = 0
    mem = 0
    for who in [resource.RUSAGE_CHILDREN, resource.RUSAGE_SELF]:
        info = resource.getrusage(who)
        utime += info.ru_utime
        stime += info.ru_stime
        mem += info.ru_maxrss
    # Memory is returned in KiB, scale to bytes
    return {"user_time": utime, "sys_time": stime, "memory": mem * 1024}


def run_benchmark_process(tool, L, queue):
    sample_size = 500
    func = tool_map[tool]
    func(sample_size=sample_size, chrom_length_mb=L)
    queue.put(process_resources())


def run_benchmark(work):

    tool, L = work
    queue = multiprocessing.Queue()
    p = multiprocessing.Process(target=run_benchmark_process, args=(tool, L, queue))
    before = time.perf_counter()
    p.start()
    perf_metrics = queue.get()
    p.join()
    if p.exitcode < 0:
        raise ValueError("Error occured", p.exitcode)
    return {"L": L, "tool": tool, **perf_metrics}


tool_map = {
    "msprime": sim_msprime,
    "hybrid": sim_msprime_hybrid,
    "ARGON": sim_argon,
}


@click.command()
@click.option("--replicates", type=int, default=5)
@click.option("--processes", type=int, default=None)
def benchmark(replicates, processes):
    """
    Runs the benchmark between msprime and ARGON.
    """

    cpu = cpuinfo.get_cpu_info()
    with open("data/dtwf_perf_cpu.txt", "w") as f:
        for k, v in cpu.items():
            print(k, "\t", v, file=f)

    chrom_lengths = np.linspace(1, 100, 20).astype(int)
    work = []
    for L in chrom_lengths:
        for name in tool_map.keys():
            work.extend([(name, L)] * replicates)

    data = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=processes) as executor:
        futures = [executor.submit(run_benchmark, item) for item in work]
        for future in concurrent.futures.as_completed(futures):
            data.append(future.result())
            print(data[-1])
            df = pd.DataFrame(data)
            df.to_csv("data/dtwf-perf.csv")


@click.command()
@click.option("--replicates", type=int, default=100)
def validate(replicates):
    """
    Validate that we are simulating the same things in the simulators
    by running some replicates and plotting the distributions of the
    number of output trees.
    """
    # NOTE: we seem to consistently get more trees from ARGON. Looking
    # at the qqplots, the distributions looks about the same, but there's
    # consistently more from ARGON. We've check the parameters as closely
    # as we can here, so I'm not sure there's much we can do.
    # However, see the discussion here:
    # https://github.com/tskit-dev/msprime-1.0-paper/pull/109
    # When we export to a tree sequence and squash the trees down properly,
    # we get the same distributions. So, this is fine.
    L = 1  # Megabases
    sample_size = 10

    nt_argon = np.zeros(replicates)
    nt_hybrid = np.zeros(replicates)
    nt_msprime = np.zeros(replicates)

    with click.progressbar(range(replicates)) as bar:
        for j in bar:
            nt_argon[j] = sim_argon(sample_size, L, count_trees=True)
            nt_hybrid[j] = sim_msprime_hybrid(sample_size, L)
            nt_msprime[j] = sim_msprime(sample_size, L)

    print(
        "mean number of trees:",
        "argon=",
        np.mean(nt_argon),
        "msprime breakpoints=",
        np.mean(nt_msprime),
        "hybrid breakpoints=",
        np.mean(nt_hybrid),
    )

    sm.graphics.qqplot(nt_argon)
    sm.qqplot_2samples(nt_argon, nt_msprime, line="45")
    plt.xlabel("argon")
    plt.ylabel("msprime")
    plt.savefig("figures/verify_argon_v_msprime.png")

    plt.close("all")


@click.group()
def cli():
    pass


cli.add_command(validate)
cli.add_command(benchmark)

if __name__ == "__main__":
    cli()
