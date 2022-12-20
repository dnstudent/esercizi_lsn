from typing import Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def plot_params_trajectory(annealing_data: pd.DataFrame, kdes: int = 2, ax: Optional[plt.Axes] = None, **kwargs):
    """Plots the trajectory (in parameter space) of the simulated annealing process.

    Args:
        annealing_data (pd.DataFrame): DataFrame of annealing steps. Must list at least the parameters and temperatures.
        kwargs: Keyworded arguments passed to DataFrame.plot.scatter.
    """
    annealing_data["logT"] = annealing_data["T"].transform("log10")
    annealing_data["catT"] = np.digitize(
        annealing_data["logT"], np.linspace(annealing_data["logT"].min(), annealing_data["logT"].max(), kdes, endpoint=False))
    annealing_data = annealing_data.rename(
        columns={"mu": r"$\mu$", "sigma": r"$\sigma$"})
    sns.scatterplot(annealing_data, x=r"$\mu$", y=r"$\sigma$",
                    hue="logT", palette="coolwarm", s=5, ax=ax)
    sns.kdeplot(annealing_data, x=r"$\mu$", y=r"$\sigma$",
                hue="catT", palette="coolwarm", legend=False, ax=ax)


def plot_energy(annealing_data: pd.DataFrame, **kwargs):
    """Plots the estimate of the energy depending on the annealing step.

    Args:
        annealing_data (pd.DataFrame): DataFrame of annealing steps, with a column "H_estimate" and "H_error".
    """
    _, ax = plt.subplots(1, 1, figsize=(15, 5))
    ax.errorbar(x="index", y="H_estimate", yerr="H_error", ecolor="red",
                data=annealing_data.iloc[1::200].reset_index(drop=False))
    return ax


def plot_histogram(data: pd.DataFrame, ax: Optional[plt.Axes] = None, **kwargs) -> plt.Axes:
    """Plots the histogram in data. Must have a column "psi" with the counts and a column "l_edge" with the left edges of the bars.

    Args:
        data (pd.DataFrame): DataFrame with the data.
        ax (Optional[plt.Axes], optional): Axes where the histogram is plotted. Defaults to None.

    Returns:
        plt.Axes: Matplotlib Axes with the plot.
    """
    if ax is None:
        _, ax = plt.subplots(1, 1)
    bin_width = data["l_edge"].iloc[1] - data["l_edge"].iloc[0]
    data = data.copy()
    #  Normalizing to compute the density histogram
    data["psi"] = data["psi"] / (data["psi"].sum()*bin_width)
    # Plotting
    ax.bar(x="l_edge", height="psi", width=bin_width,
           align='edge', data=data, **kwargs)
    return ax
