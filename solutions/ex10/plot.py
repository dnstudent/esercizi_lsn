from pathlib import Path
import matplotlib.pyplot as plt
from typing import List


import cartopy.crs as ccrs
import cartopy.feature as cfeature

from global_utils import plot_trajectory
from .data import *
from .paths import *


def plot_genetic_trajectory(algo: Algo, coordinates: Path, p_line: int,  ax: plt.Axes):
    """Plots the trajectory of the best result in a TSP.

    Args:
        tsp (TSP): TSP which was run.
        algo (Algo): Algorithm which was used for crossover.
        coordinates (Path): Path to the cities coordinates.
        p_line (int): Line number in the primes file.
        ax (plt.Axes): Matplotlib axes where the result will be plotted.
    """
    results = data_from(algo, p_line)
    coords = pd.read_csv(coordinates, names=["longitude", "latitude"])
    sorted_results = results.sort_values("total_distance")
    best_fitness = sorted_results["total_distance"].iloc[0]
    best_individual = sorted_results.drop(columns="total_distance").iloc[0]
    coords_sequence = coords.loc[best_individual]
    coords_sequence["t"] = range(len(coords_sequence))
    plot_trajectory("longitude", "latitude", ax,
                    data=coords_sequence, cmap="plasma")
    ax.scatter("longitude", "latitude", s=4, data=coords)
    ax.set_title(f"Total distance: {best_fitness:.3f}")


def plot_panel(input: Path, seeds: List[int]):
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
            plot_genetic_trajectory(algo, input, p_line=seed, ax=ax)
            ax.set_xticklabels([])
            ax.set_yticklabels([])
    for ax, seed in zip(axes, seeds):
        ax.set_ylabel(f"seed: {seed}")
        ax.yaxis.set_label_position("right")
    return fig


def plot_capital_names(item, ax):
    ax.text(
        item["longitude"], item["latitude"], item["city"], transform=ccrs.Geodetic()
    )
    return ax


def plot_route(individual, capitals, ax):
    route = capitals.loc[individual]
    start = route.iloc[0]
    finish = route.iloc[-1]
    ax.plot("longitude", "latitude", data=route, transform=ccrs.Geodetic())
    ax.plot(
        "longitude",
        "latitude",
        c="r",
        marker="o",
        data=start,
        transform=ccrs.Geodetic(),
    )
    ax.text(start["longitude"], start["latitude"],
            "Start", transform=ccrs.Geodetic())
    ax.plot(
        "longitude",
        "latitude",
        c="g",
        marker="o",
        data=finish,
        transform=ccrs.Geodetic(),
    )
    ax.text(
        finish["longitude"], finish["latitude"], "Finish", transform=ccrs.Geodetic()
    )
    return ax

def plot_on_map(algo, n_continents, capitals, fusion_p: Optional[float]):
    fig = plt.figure(dpi=300, figsize=(18, 9))
    ax = fig.add_subplot(projection=ccrs.LambertConformal())
    ax.set_extent([-160, -65, 22, 60])
    ax.gridlines(draw_labels=True, dms=True, x_inline=False,
                y_inline=False, linestyle="--")
    ax.add_feature(cfeature.STATES)
    ax.scatter(
        "longitude", "latitude", marker="*", s=10, data=capitals, transform=ccrs.Geodetic()
    )

    best_route = data_from(algo, n_continents, fusion_p).sort_values("total_distance").drop(columns="total_distance").iloc[0]
    return plot_route(best_route, capitals, ax)