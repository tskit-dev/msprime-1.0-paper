import os
import sys

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

panels = [["coalescent", "dtwf"],
          ["smc", "sweep"]]

palette = {
           "ms": "#F7C362",
           "msprime": "#2E8CB9",
           "msms": "#8DA470",
           "discoal": "#30908E",
           "argon": "#Df4D47",
           "msprime_smc": "#66A9FF",
           "scrm": "#9E5De8",
          }

params = {'xtick.labelsize': 9,
          'ytick.labelsize': 9}

def plot_all(data_dir):
    sns.set(style="darkgrid")
    plt.rcParams.update(params)

    rows = len(panels)
    cols = len(panels[0])
    fig, axes = plt.subplots(rows, cols, squeeze=False)

    for row in range(rows):
        for col in range(cols):
            sim = panels[row][col]
            data = load(f"{data_dir}/{sim}")

            ax = axes[row][col]
            sns.lineplot(x="sample_size", y="time", hue="tool", palette=palette,
                         data=data, ax=ax)
            remove_subplot_legend_title(ax);
            ax.set_title(sim)
            ax.set_xlabel("Sample Size" if row > 0 else "")
            ax.set_ylabel("Time (s)" if col == 0 else "")

    save(fig)

def remove_subplot_legend_title(ax):
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles=handles[1:], labels=labels[1:])

def load(sim_path):
    return pd.read_csv(f"{sim_path}.csv")

def save(fig):
    plt.tight_layout()
    fig.savefig(get_output_filename())

def get_output_filename():
    basedir = "plots"
    if not os.path.exists(basedir):
        os.mkdir(basedir)
    return os.path.join(basedir, "plot.png")

if __name__ == "__main__":
    plot_all(sys.argv[1])
