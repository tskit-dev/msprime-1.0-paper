import click
import tskit
import msprime


@click.command()
def mutated_tree():
    """
    Make a figure with (a) a tree and (b) some mutations added to it.
    """
    ts = msprime.sim_ancestry(
        50,
        population_size=1e4,
        recombination_rate=0.5e-8,
        sequence_length=1000,
        random_seed=14,
    )

    model = msprime.F84(kappa=2)
    mts = msprime.sim_mutations(ts, rate=5e-8, model=model, random_seed=5)


    def do_svg(ts, **kwargs):
        ts.draw_svg(
            size=(300, 200),
            node_labels={},
            mutation_labels={m.id: m.derived_state for m in ts.mutations()},
            symbol_size=5,
            force_root_branch=True,
            **kwargs
        )


    do_svg(ts, path="illustrations/unmutated_tree.svg")
    do_svg(mts, path="illustrations/mutated_tree.svg")


@click.command()
def arg_ts():
    tables = tskit.TableCollection(1.0)
    tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, time=0.5)
    # tables.nodes.add_row(flags=msprime.NODE_IS_RE_EVENT, time=0.5)
    tables.nodes.add_row(flags=0, time=1.0)
    tables.nodes.add_row(flags=0, time=1.5)
    tables.nodes.add_row(flags=0, time=2.0)

    tables.edges.add_row(0, 1, 5, 0)
    tables.edges.add_row(0, 1, 3, 1)
    tables.edges.add_row(0.0, 0.3, 5, 3)
    tables.edges.add_row(0.3, 1, 4, 3)
    tables.edges.add_row(0, 1, 4, 2)
    tables.edges.add_row(0, 1, 6, 5)
    tables.edges.add_row(0, 1, 6, 4)

    tables.sort()
    print(tables)
    ts = tables.tree_sequence()

    # https://stackoverflow.com/questions/46077392/additional-options-in-chrome-headless-print-to-pdf
    # Sizes are probably all wrong here
    style = """\
        @media print {
          @page { margin: 0; size: 4in 3in}
          body { margin: 1.6cm; }
        }
    """
    svgfile = "illustrations/arg-ts.svg"
    ts.draw_svg(svgfile, size=(400, 300), style=style)

    print(ts.draw_text())

    ts = ts.simplify()
    svgfile = "illustrations/arg-ts-simplified.svg"
    ts.draw_svg(svgfile, size=(400, 300), style=style)
    print(ts.draw_text())


@click.group()
def cli():
    pass


cli.add_command(mutated_tree)
cli.add_command(arg_ts)

if __name__ == "__main__":
    cli()
