import subprocess
import sys
import time

import msprime

REPLICATES = 1

Ne = 1e4
mu = 1e-9
r = 1e-8

def run_all():
    run_standard_coalescent()
    run_dtwf()
    run_smc()
    run_sweep()


def run_standard_coalescent():
    L = mb(1)
    sample_sizes = range(1000, 10000 + 1, 1000)

    def ms(sample_size):
        return ms_run(sample_size, L)
    def msms(sample_size):
        return msms_run(sample_size, L)
    def msp(sample_size):
        return msp_run(sample_size=sample_size, length=L)

    create_data(sample_sizes,
                [ms, msms, msp],
                ["ms", "msms", "msprime"],
                "coalescent")


def run_dtwf():
    L = mb(1)
    sample_sizes = range(1000, 10000 + 1, 1000)

    def argon(sample_size):
        return argon_run(sample_size, Ne, L)
    def msp(sample_size):
        return msp_run(sample_size=sample_size, length=L, Ne=Ne, model="dtwf")

    create_data(sample_sizes, [msp, argon], ["msprime", "argon"], "dtwf")


def run_smc():
    L = mb(10)

    sample_sizes = range(10000, 20000 + 1, 1000)
    def scrm(sample_size):
        return scrm_run(sample_size, L)
    def msprime_exact(sample_size):
        return msp_run(sample_size=sample_size, length=L)
    def msprime_smc(sample_size):
        return msp_run(sample_size=sample_size, length=L, model="smc")

    create_data(sample_sizes,
                [scrm, msprime_exact, msprime_smc],
                ["scrm", "msprime", "msprime_smc"],
                "smc")


# TODO Set up sweep parameters.
def run_sweep():
    L = mb(1)

    sample_sizes = range(1000, 15000 + 1, 1000)

    def msms(sample_size):
        return msms_run(sample_size, L)
    def discoal(sample_size):
        return discoal_run(sample_size, L)
    def msp(sample_size):
        return msp_run(sample_size=sample_size, length=L)

    create_data(sample_sizes,
                [msms, discoal, msp],
                ["msms", "discoal", "msprime"],
                "sweep")


def ms_run(sample_size, L):
    args = ms_style_args(sample_size, L)
    exec_cli("./sims/ms", args)


def argon_run(sample_size, Ne, L):
    length_mb = L / mb(1)
    # ARGON assumes a haploid Ne whereas msprime uses diploid
    haploid_ne = 2 * int(Ne)
    args = (
        f"-N {haploid_ne} -pop 1 {sample_size} "
        f"-size {length_mb} -rec {r} -quiet -screen"
    )

    for _ in range(REPLICATES):
        exec_jar("sims/ARGON.jar", args)


def scrm_run(sample_size, L):
    ms_args = ms_style_args(sample_size, L)
    scrm_args = "-l 0"
    args = ms_args + " " + scrm_args
    exec_cli("./sims/scrm", args)


def msms_run(sample_size, L):
    ms_args = ms_style_args(sample_size, L)
    args = f"-N {Ne} -ms {ms_args}"
    exec_jar("sims/msms.jar", args)


def discoal_run(sample_size, L):
    ms_args = ms_style_args(sample_size, L)
    args = ms_args + " " + f"{int(L)}"
    exec_cli("./sims/discoal", args)


def msp_run(**kwargs):
    if "Ne" not in kwargs:
        kwargs["Ne"] = Ne
    if "recombination_rate" not in kwargs:
        kwargs["recombination_rate"] = r
    kwargs["num_replicates"] = REPLICATES

    replicates = msprime.simulate(**kwargs)
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
    with open(f"data/{fname}.csv", "w") as f:
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
    return (end - start) / REPLICATES


def mb(L):
    return L * 1e6


if __name__ == "__main__":
    if "hudson" in sys.argv:
        run_standard_coalescent()
    if "dtwf" in sys.argv:
        run_dtwf()
    if "smc" in sys.argv:
        run_smc()
    if "sweep" in sys.argv:
        run_sweep()

    if len(sys.argv) == 1:
        run_all()
