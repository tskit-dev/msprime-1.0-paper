import click
import tskit
import msprime

import matplotlib.pyplot as plt


@click.command()
def mutated_tree():
    """
    Make a figure with (a) a tree and (b) some mutations added to it.
    """
    ts = msprime.sim_ancestry(
        3,
        population_size=1e4,
        recombination_rate=1e-8,
        sequence_length=1000,
        random_seed=96,
    )

    model = msprime.F84(kappa=2)
    mts = msprime.sim_mutations(ts, rate=1e-7, model=model, random_seed=4)

    height = 280
    width = 700
    top = 50

    def do_svg(ts, **kwargs):
        colours = plt.rcParams['axes.prop_cycle'].by_key()['color']

        # Sizes are probably all wrong here
        style = """\
            @media print {
              @page { margin: 0; size: 6in 2.5in}
              body { margin: 1.6cm; }
            }
            """
        for j in range(ts.num_individuals):
            style += f"\n.node.i{j} > .sym {{fill: {colours[j]}}}"

        mut_colour = colours[ts.num_individuals]
        style += (
            f".mut text {{fill: {mut_colour}; font-style: italic}}"
            f".mut .sym {{fill: none; stroke: {mut_colour}}}"
        )

        return ts.draw_svg(
            size=(width / 2, height - top),
            node_labels={},
            mutation_labels={m.id: m.derived_state for m in ts.mutations()},
            symbol_size=5,
            style=style,
            **kwargs,
        )

    font_size = 15

    # I think serif is the default, and matches what we're using for labels?
    def make_text(text, y, font_family="serif"):
        return (
            f'<text x="{width / 4}" y="{y}" font-size="{font_size}" '
            f'font-family="{font_family}" text-anchor="middle">'
            f"{text}</text>"
        )

    svg1 = do_svg(ts)
    svg2 = do_svg(mts)
    fig = (
        f'<svg baseProfile="full" height="{height+top}" version="1.1" width="{width}" '
        'xmlns="http://www.w3.org/2000/svg" xmlns:ev="http://www.w3.org/2001/xml-events" '
        'xmlns:xlink="http://www.w3.org/1999/xlink">'
        f'<g transform="translate(0 {top})">'
        + make_text("(A)", y=-20)
        + make_text("sim_ancestry", y=-5, font_family="monospace")
        + svg1
        + "</g>"
        f'<g transform="translate({width/2} {top})">'
        + make_text("(B)", y=-20)
        + make_text("sim_mutations", y=-5, font_family="monospace")
        + svg2
        + "</g>"
        "</svg>"
    )
    with open("illustrations/mutated_tree.svg", "w") as f:
        f.write(fig)


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
