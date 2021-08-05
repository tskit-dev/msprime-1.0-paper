"""
Generate data for the ARG figure.
"""
import numpy as np
import msprime
import pandas as pd


def run_arg_sim():

    L_col = []
    size_col = []
    arg_col = []
    for megabases in np.linspace(0.1, 5, 20):
        L = int(megabases * 1_000_000)
        arg_ts = msprime.sim_ancestry(
            100,
            population_size=10_000,
            sequence_length=L,
            recombination_rate=1e-8,
            random_seed=42,
            record_full_arg=True,
        )
        flags = arg_ts.tables.nodes.flags
        # Samples have flags == 1 and ordinary coalescent nodes have
        # flags == 0. So, anything > 1 is an ARG node.
        arg_nodes = flags > 1
        L_col.append(L)
        arg_fraction = np.sum(arg_nodes) / arg_ts.num_nodes
        ts = arg_ts.simplify()
        size_ratio = ts.tables.nbytes / arg_ts.nbytes
        size_col.append(size_ratio)
        arg_col.append(arg_fraction)
        print(L, arg_fraction, size_ratio)
    data = {"L": L_col, "size_ratio": size_col, "arg_nodes": arg_col}

    df = pd.DataFrame(data)
    print(df)
    df.to_csv("data/arg.csv")


if __name__ == "__main__":
    run_arg_sim()
