from pathlib import Path
import matplotlib.pyplot as plt
from typing import List

from global_utils import plot_trajectory
from .data import *
from .paths import *


def plot_genetic_trajectory(tsp: TSP, algo: Algo, coordinates: Path, p_line: int,  ax: plt.Axes):
    """Plots the trajectory of the best result in a TSP.

    Args:
        tsp (TSP): TSP which was run.
        algo (Algo): Algorithm which was used for crossover.
        coordinates (Path): Path to the cities coordinates.
        p_line (int): Line number in the primes file.
        ax (plt.Axes): Matplotlib axes where the result will be plotted.
    """
    results = data_from(tsp, algo, p_line)
    coords = pd.read_csv(coordinates, names=["x", "y"])
    sorted_results = results.sort_values("total_distance")
    best_fitness = sorted_results["total_distance"].iloc[0]
    best_individual = sorted_results.drop(columns="total_distance").iloc[0]
    coords_sequence = coords.loc[best_individual]
    coords_sequence["t"] = range(len(coords_sequence))
    plot_trajectory("x", "y", ax, data=coords_sequence, cmap="plasma")
    ax.scatter("x", "y", s=4, data=coords)
    ax.set_title(f"Total distance: {best_fitness:.3f}")


def plot_panel(tsp: TSP, input: Path, seeds: List[int]):
    """Plots a panel to compare the results given by different algorithms and different seeds.

    Args:
        tsp (TSP): TSP which was run.
        input (Path): Path to the cities coordinates.
        seeds (List[int]): List of seeds which was used in the runs.
    """
    # Preparing the matplotlib figure
    fig = plt.figure(figsize=(len(Algo)*3, len(seeds)*3),
                     constrained_layout=True)
    # Making n_algos subfigures (the panel's columns)
    subfigs = fig.subfigures(1, len(Algo))
    fig.suptitle("Best results from each algorithm and seed")

    for algo, subfig in zip(Algo, subfigs):
        #  Giving each subfigures the algorithm's name as title
        subfig.suptitle(f"{algo.value}")
        axes = subfig.subplots(len(seeds), 1, sharex=True, sharey=True,
                               subplot_kw={"aspect": "equal"})
        for seed, ax in zip(seeds, axes):
            plot_genetic_trajectory(tsp, algo, input, p_line=seed, ax=ax)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
    for ax, seed in zip(axes, seeds):
        ax.set_ylabel(f"seed: {seed}")
        ax.yaxis.set_label_position("right")
    return fig
