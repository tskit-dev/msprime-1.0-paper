import subprocess
import sys
import time


#import sys
#sys.path.insert(0, "/home/franz/syncedprojects/msprime_performance/msprime")

import msprime

REPLICATES = 1

Ne = 1
mu = 1e-9
r = 0


def run_simbac():
    L = mb(4.5)
    #L = 10000
    #sample_sizes = list(range(10, 100 + 1, 10)) + list(range(200, 900 + 1, 100)) + list(range(1000, 10000 + 1, 1000))
    sample_sizes = list(range(10,10+1,1)) + list(range(50,500 + 1, 50))
    #sample_sizes = range(10, 100 + 1, 10)
    #sample_sizes = range(200, 900 + 1, 100)

    def simbac(sample_size):
        return simbac_run(sample_size, L)
    def msp(sample_size):
        return msp_run(samples=sample_size, sequence_length=L, population_size=1,
                       ploidy = 1, recombination_rate=0,
                       gene_conversion_rate=0.015, gene_conversion_tract_length=500)

    create_data(sample_sizes, [msp, simbac], ["msprime", "simbac"], "simbac")


def ms_run(sample_size, L):
    args = ms_style_args(sample_size, L)
    exec_cli("./sims/ms", args)


def simbac_run(sample_size, L):
    ms_args = f"-N {sample_size} -B {int(L)}"
    # using R=2*0.015 as gene conversion/recombination rate as SimBac uses R/2
    simbac_args = "-R 0.03"
    args = ms_args + " " + simbac_args
    for _ in range(REPLICATES):
        exec_cli("./sims/SimBac", args)


def msp_run(**kwargs):
    if "population_size" not in kwargs:
        kwargs["population_size"] = Ne
    if "recombination_rate" not in kwargs:
        kwargs["recombination_rate"] = r
    kwargs["num_replicates"] = REPLICATES

    replicates = msprime.sim_ancestry(**kwargs)
    for ts in replicates:
        pass


def exec_cli(command, args_str):
    arg_list = [command] + args_str.split()
    print(" ".join(arg_list))
    subprocess.call(arg_list, stdout=subprocess.DEVNULL)


def exec_jar(jar_path, sim_args_str):
    args = f"-jar {jar_path} {sim_args_str}"
    exec_cli("java", args)


def ms_style_args(sample_size, L):
    rho = 4 * Ne * r * L
    return f"{sample_size} {REPLICATES} -t 10 -r {rho} {int(L)}"


def create_data(sample_sizes, tools, tool_names, fname):
    runtimes = run_sims(sample_sizes, tools, tool_names)
    with open(f"{fname}.csv", "w") as f:
        f.write(f"sample_size,tool,time\n")
        for sample_size, runtime_row in zip(sample_sizes, runtimes):
            for tool_name, tool_time in zip(tool_names, runtime_row):
                row = f"{sample_size},{tool_name},{tool_time}\n"
                f.write(row)
        f.close()


def run_sims(sample_sizes, tools, tool_names):
    runtimes = []
    for sample_size in sample_sizes:
        times = []
        for tool, name in zip(tools, tool_names):
            print(f"Running {name} on {sample_size} samples...")
            times.append(time_tool(tool, sample_size))
        runtimes.append(times)
    return runtimes


def time_tool(tool, sample_size):
    start = time.time()
    tool(sample_size)
    end = time.time()
    print(f"Took {round(end - start, 3)} seconds to complete {REPLICATES} replicates")
    return (end - start) / REPLICATES


def mb(L):
    return L * 1e6


if __name__ == "__main__":
    if "simbac" in sys.argv:
        run_simbac()
