import matplotlib
import numpy as np
import pandas as pd
import msprime
import attr
import time
import sys
import functools
import subprocess
import matplotlib
from matplotlib import pyplot as plt

"""
this outputs a csv file of the performance comparison (time_comp.csv)
along with a rough figure to visualize stuff (time_comp.pdf)
"""

_mspms_executable = [sys.executable, "mspms_dev.py"]
_discoal_executable = ["../tools/discoal"]


class Timer:
    def __init__(self):
        self._start_time = None

    def start(self):
        if self._start_time is not None:
            raise Exception(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        if self._start_time is None:
            raise Exception(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None
        return(elapsed_time)

@attr.s
class Run:
    """
    The superclass of all tests. The only attribute defined is the output
    directory for the test, which is guaranteed to exist when the
    test method is called.
    """

    output_dir = attr.ib(type=str, default=None)

    def _run_simulation(self, args):
        p1 = subprocess.Popen(args, stdout=subprocess.PIPE)
        output = p1.communicate()[0]
        p1.wait()
        if p1.returncode != 0:
            raise ValueError("Error occured in subprocess: ", p1.returncode)

    def get_ms_seeds(self):
        max_seed = 2 ** 16
        seeds = [np.random.randint(1, max_seed) for j in range(3)]
        return ["-seed"] + list(map(str, seeds))

    def _run_msprime(self, args):
        return self._run_simulation(
            _mspms_executable + args.split() + self.get_ms_seeds()
        )

class DiscoalRun(Run):
    def get_discoal_seeds(self):
        max_seed = 2 ** 16
        seeds = [np.random.randint(1, max_seed) for j in range(3)]
        return ["-d"] + list(map(str, seeds))

    def _discoal_str_to_ms(self, args):
        # convert discoal string to msprime string
        tokens = args.split(" ")
        # cut out sites param
        del tokens[2]
        # adjust popIDs
        for i in range(len(tokens)):
            # pop size change case
            if tokens[i] == "-en":
                tokens[i + 2] = str(int(tokens[i + 2]) + 1)
            # migration rate case
            if tokens[i] == "-m":
                tokens[i + 1] = str(int(tokens[i + 1]) + 1)
                tokens[i + 2] = str(int(tokens[i + 2]) + 1)
        msp_str = " ".join(tokens)
        return msp_str

    def _run_discoal(self, args):
        return self._run_simulation(
            _discoal_executable + args.split() + self.get_discoal_seeds()
        )


    def _discoal_str_to_simulation(self, args):
        # takes discoal command line as input
        # and returns msprime run treeseqs

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
                raise ValueError(
                    "sweeps with population size changes remain unimplemented"
                )
            # migration rate case
            if (tokens[i] == "-m") or (tokens[i] == "-p"):
                raise ValueError(
                    "sweeps with multiple populations remain unimplemented"
                )
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
        replicates = msprime.sim_ancestry(
            [msprime.SampleSet(sample_size, ploidy=1)],
            population_size=refsize,
            model=mod_list,
            recombination_rate=recomb_rate,
            sequence_length=seq_length,
            discrete_genome=False,
            num_replicates=nreps,
        )
        mutate = functools.partial(
            msprime.sim_mutations, discrete_genome=False, rate=mu
        )
        return map(mutate, replicates)


class DiscoalSweeps(DiscoalRun):
    """
    Compare the result of sweeps in msprime and discoal.
    """

    def _run(self, args):
        t = Timer()
        t.start()
        replicates = self._discoal_str_to_simulation(args)
        time1 = t.stop()
        t.start()
        self._run_discoal(args)
        time2 = t.stop()
        return time1, time2

def main():
    ds = DiscoalSweeps()
    msp = []
    dsc = []
    r = 1e-9
    u = 1e-9
    refsize = 1e6
    seqlen = np.linspace(200,200000, 10) 
    for l in seqlen:
        rho = 4 * refsize * (l-1) * r
        theta = 4 * refsize * (l-1) * u
        cmd = f"10 100 {int(l)} -t {theta} -r {rho} -ws 0 -a 1000 -x 0.5 -N 10000"
        tup = ds._run(cmd)
        msp.append(tup[0])
        dsc.append(tup[1])
    df = pd.DataFrame({'msprime':msp,
                       'discoal':dsc,
                       'sequence length':seqlen
                      })
    df.set_index('sequence length', inplace=True)
    df.to_csv("time_comp.csv")
    fig = df.plot().get_figure()
    fig.savefig('time_comp.pdf')

if __name__ == "__main__":
    main()
