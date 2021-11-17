import subprocess
import concurrent.futures
import multiprocessing
import resource

import numpy as np
import click
import msprime
import pandas as pd
import cpuinfo


def _discoal_str_to_msprime(args):
    # takes discoal command line as input and returns an iterator over the
    # msprime tree sequences.

    tokens = args.split(" ")
    # positional args
    sample_size = int(tokens[0])
    nreps = int(tokens[1])
    seq_length = int(tokens[2])
    # parse discoal command line for params
    # init ones we definitely need for comparison
    theta = rho = alpha = sweep_site = sweep_mod_time = None
    refsize = 1e6
    for i in range(3, len(tokens)):
        # pop size change case
        if tokens[i] == "-en":
            raise ValueError("sweeps with population size changes remain unimplemented")
        # migration rate case
        if (tokens[i] == "-m") or (tokens[i] == "-p"):
            raise ValueError("sweeps with multiple populations remain unimplemented")
        # split or admixture case
        if (tokens[i] == "-ea") or (tokens[i] == "-ed"):
            raise ValueError("sweeps with splits or admixture not supported")
        # sweep params
        if tokens[i] == "-x":
            sweep_site = float(tokens[i + 1])
        if (tokens[i] == "-ws") or (tokens[i] == "-wd") or (tokens[i] == "-wn"):
            sweep_mod_time = float(tokens[i + 1])
        if tokens[i] == "-a":
            alpha = float(tokens[i + 1])
        if tokens[i] == "-N":
            refsize = float(tokens[i + 1])
        # coalescent params
        if tokens[i] == "-t":
            theta = float(tokens[i + 1])
        if tokens[i] == "-r":
            rho = float(tokens[i + 1])
    mod_list = []
    if alpha is not None:
        # sweep model
        s = alpha / (2 * refsize)
        mod = msprime.SweepGenicSelection(
            position=np.floor(sweep_site * seq_length),
            start_frequency=1.0 / (2 * refsize),
            end_frequency=1.0 - (1.0 / (2 * refsize)),
            s=s * 2,  # discoal fitness model is 1, 1+s, 1+2s
            dt=1e-6,
        )
        mod_list.append(msprime.StandardCoalescent(duration=sweep_mod_time))
        mod_list.append(mod)
        # if an event is defined from discoal line
        # best thing to do is rescale to Ne=0.25
        # so that time scale are consistent
        # see note at msprime/cli.py line 626
        # and following for alternate solution
        if sweep_mod_time > 0:
            refsize = 0.25
            mod.s = alpha / refsize
    # append final model
    mod_list.append("hudson")
    # scale theta and rho
    recomb_rate = rho / (4 * refsize * (seq_length - 1))
    mu = theta / (4 * refsize * seq_length)
    # We're only interested in ancestry sim here.
    assert mu == 0
    replicates = msprime.sim_ancestry(
        [msprime.SampleSet(sample_size, ploidy=1)],
        population_size=refsize,
        model=mod_list,
        recombination_rate=recomb_rate,
        sequence_length=seq_length,
        discrete_genome=False,
        num_replicates=nreps,
    )
    return replicates


def run_discoal(args):
    cmd = "./tools/discoal " + args
    subprocess.run(cmd, check=True, shell=True, stdout=subprocess.DEVNULL)


def run_msprime(args):
    replicates = _discoal_str_to_msprime(args)
    for _ in replicates:
        pass

def get_process_resources():
    utime = 0
    stime = 0
    mem = 0
    for who in [resource.RUSAGE_CHILDREN, resource.RUSAGE_SELF]:
        info = resource.getrusage(who)
        utime += info.ru_utime
        stime += info.ru_stime
        mem += info.ru_maxrss
    # Memory is returned in KiB on Linux, scale to bytes
    return {"user_time": utime, "sys_time": stime, "memory": mem * 1024}


def run_benchmark_process(tool, cmd, queue):
    func = tool_map[tool]
    func(cmd)
    queue.put(get_process_resources())


def run_benchmark(work):

    tool, L, cmd = work
    queue = multiprocessing.Queue()
    p = multiprocessing.Process(target=run_benchmark_process, args=(tool, cmd, queue))
    p.start()
    perf_metrics = queue.get()
    p.join()
    if p.exitcode < 0:
        raise ValueError("Error occured", p.exitcode)
    return {"L": L, "tool": tool, **perf_metrics}


tool_map = {
    "msprime": run_msprime,
    "discoal": run_discoal,
}


@click.command()
@click.option("--replicates", type=int, default=3)
@click.option("--processes", type=int, default=1)
def benchmark(replicates, processes):
    """
    Runs the benchmark between msprime and discoal
    """

    cpu = cpuinfo.get_cpu_info()
    with open("data/sweeps_perf_cpu.txt", "w") as f:
        for k, v in cpu.items():
            print(k, "\t", v, file=f)

    r = 1e-8
    # Set theta to 0 to avoid conflating mutation generation
    theta = 0
    refsize = 1e4
    # 100kb up to 3000kb
    Ls = np.linspace(100_000, 3_000_000, 20).astype(int)
    work = []
    for L in Ls:
        for name in tool_map.keys():
            rho = 4 * refsize * (L - 1) * r
            cmd = (
                f"200 {replicates} {L} -t {theta} -r {rho} -ws 0 -a 50 -x 0.5 -N {refsize}"
            )
            work.extend([(name, L, cmd)])

    data = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=processes) as executor:
        futures = [executor.submit(run_benchmark, item) for item in work]
        for future in concurrent.futures.as_completed(futures):
            data.append(future.result())
            print(data[-1])
            df = pd.DataFrame(data)
            df['user_time'] = df.user_time / replicates
            df['sys_time'] = df.sys_time / replicates
            df.to_csv("data/sweeps_perf.csv")


@click.group()
def cli():
    pass


cli.add_command(benchmark)

if __name__ == "__main__":
    cli()
