import time

import msprime
import numpy as np

def run(L, N, num_samples):
    start = time.process_time()
    ts = msprime.sim_ancestry(
            samples=num_samples,
            population_size=N,
            sequence_length=L,
            recombination_rate=1e-8)
    end = time.process_time()
    return (max(ts.tables.nodes.time), ts.num_trees, end - start)


def csv(x):
    return ",".join(map(str, x)) + "\n"


def run_sims():
    outfile = "data/ancestry-perf.csv"
    num_reps = 3
    num_samples = [1000, 100000]
    L_VALS = [1e6, 5e6, 1e7, 5e7, 1e8]
    N_VALS = [1000, 5000, 10000, 50000, 100000, 200000, 300000]
    ln_pairs = [(L, N) for L in L_VALS for N in N_VALS]
    ln_pairs.sort(key=lambda x: x[0] * x[1])

    with open(outfile, "w") as f:
        f.write(csv(["N", "L", "num_samples", "height", "num_trees", "time"]))
        for L, N in ln_pairs:
            for ns in num_samples:
                if L * N * np.log10(ns) < 3e12:
                    for _ in range(num_reps):
                        print(f"... L={L}, N={N}, num_samples={ns}, L*N={int(L*N)}.")
                        h, n, t = run(L, N, ns)
                        f.write(csv([N, L, ns, h, n, t]))
                        f.flush()

if __name__ == "__main__":
    run_sims()
