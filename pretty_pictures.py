#!/usr/bin/env python

import msprime


## Adding mutations
# 
# The goal here is to make a figure with (a) a tree and (b) some mutations added to it.


ts = msprime.sim_ancestry(
        50,
        population_size=1e4,
        recombination_rate=0.5e-8,
        sequence_length=1000,
        random_seed=14
)

model = msprime.F84(kappa=2)
mts = msprime.sim_mutations(
        ts,
        rate=5e-8,
        model=model,
        random_seed=5
)


def do_svg(ts, **kwargs):
    ts.draw_svg(
            size=(300,200),
            node_labels={},
            mutation_labels={m.id: m.derived_state for m in ts.mutations()},
            symbol_size=5,
            force_root_branch=True,
            **kwargs
    )

do_svg(ts, path="figures/unmutated_tree.svg")
do_svg(mts, path="figures/mutated_tree.svg")

