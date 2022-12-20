from typing import Tuple, Optional, List, Any

import matplotlib.pyplot as plt

from .vars import *
from .data import *


# used for display pourposes
RENAMES = {"H": r"$\mathrm{H}$", "Sum_s": r"$\sum s_i$", "Sum_s2": r"$\sum s_is_j$"}


def show_equilibrations(
    T: float,
    h: float,
    samplers: Union[Sampler, List[Sampler]],
    take: Optional[int] = None,
    **kwargs,
):
    """Plots the results of an equilibration run

    Args:
        T (float): Temperature
        h (float): External magnetic field. Used only when computing m.
        samplers (Union[Sampler, List[Sampler]]): Samplers to use.
        take (Optional[List[int]], optional): How many steps to take for each sampler. If None take all of them. Defaults to None.

    Returns:
        _type_: _description_
    """
    if not hasattr(samplers, "__iter__"):
        samplers = [samplers]

    fig, axes = plt.subplots(
        3,
        len(samplers),
        figsize=(5 * len(samplers), 5),
        constrained_layout=True,
        sharex="col",
        sharey="row",
        **kwargs,
    )
    fig.suptitle(f"Equilibration, T={T}")

    dfs = read_equilibrations(T, h, sampler=samplers, nrows=take)

    def _plot(sampler, axes):
        for i, df0 in enumerate(dfs[sampler]):
            df0 = df0.rename(columns=RENAMES)
            for var, ax in zip(df0.columns, axes):
                ax.plot(
                    var,
                    data=df0,
                    label=f"{i}",
                )
                ax.set_title(f"{var} {sampler.value}")
                ax.grid(visible=True, which="major", linestyle="-", linewidth=1.5)
                ax.grid(visible=True, which="minor", linestyle="--", linewidth=1)
            ax.legend()

    for i, sampler in enumerate(samplers):
        _plot(sampler, axes[:, i])

    return axes


def show_autocorrelations(
    T: float,
    h: float,
    samplers: Union[Sampler, List[Sampler]],
    lags: Optional[int] = None,
) -> List[plt.Axes]:
    """Shows the autocorrelation results of a previous equilibration run.

    Args:
        T (float): Temperature
        h (float): External magnetic field. Used only when computing m.
        samplers (Union[Sampler, List[Sampler]]): Samplers to use.
        lags (Optional[int], optional): Number of lags for which autocorrelation will be computed. If None use every lag from 0 to n_measures. Defaults to None.
    """
    if not hasattr(samplers, "__iter__"):
        samplers = [samplers]

    dfs = read_autocorrelation(T, h, sampler=samplers)
    _, axes = plt.subplots(
        1,
        len(samplers),
        figsize=(5 * len(samplers), 5),
        squeeze=False,
        constrained_layout=True,
    )
    for (sampler, df), ax in zip(dfs.items(), axes.flatten()):
        if lags is None:
            lags = len(df)
        df.iloc[:lags].rename(columns=RENAMES).plot(ax=ax)
        ax.axhline(np.exp(-2), c="red", label="$\pm e^{-2}$")
        ax.axhline(-np.exp(-2), c="red")
        ax.set_title(f"{sampler.value} autocorrelation")
        ax.legend()
    return axes


def show_measures(
    J: float,
    T: float,
    h: float,
    N: int,
    samplers: Union[Sampler, List[Sampler]],
    **kwargs,
) -> np.ndarray[Any, plt.Axes]:
    """Plots the variables estimations of a system at temperature T.
    Must be used after "run(T, ...)".

    Args:
        J (float): Spin coupling.
        T (float): Temperature
        h (float): External magnetic field. Used only when computing m.
        N (int): Number of spins.
        samplers (Union[Sampler, List[Sampler]]): The MCMC algorithm(s) used to sample points.
        **kwargs: keyworded arguments passed to plt.errorbars()
    Returns:
        np.ndarray[Any, plt.Axes]: Array of the matplotlib axes on which data is plotted.
    """
    if type(samplers) is Sampler:
        samplers = [samplers]
    estimations = read_measures(T, h, sampler=samplers)
    theory = compute_theory(J, T, h, N)
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for variable, ax in zip(VARIABLES, axes.flatten()):
        for sampler, df in estimations.items():
            ax.errorbar(
                data=df[variable].reset_index(drop=False),
                x="index",
                y="estimate",
                yerr="error",
                label=f"{sampler.value}",
                **kwargs,
            )
        ax.axhline(theory[variable], c="red", label="theory")
        ax.set_title(variable)
    ax.legend()
    return axes


def show_results(
    J: float,
    Ts: List[float],
    h: float,
    n_spins: int,
    samplers: Union[Sampler, List[Sampler]],
) -> np.ndarray[Any, plt.Axes]:
    """Plots the difference between estimations and theory for a group of temperatures.

    Args:
        J (float): Spin coupling.
        Ts (List[float]): Temperatures.
        h (float): External magnetic field. Applied only when estimating m.
        n_spins (int): Number of spins.
        samplers (Union[Sampler, List[Sampler]]): MCMC sampling algorithm.

    Returns:
        np.ndarray[Any, plt.Axes]: The figure's axes.
    """
    results = read_results(Ts, h, samplers)
    fig, axes = plt.subplots(
        2, 2, constrained_layout=True, sharex=True, figsize=(10, 7)
    )
    fig.suptitle("Results: estimate - theory")
    theory = compute_theory(J, Ts, h, n_spins)
    for var, ax in zip(VARIABLES, axes.flatten()):
        for sampler, df in results.items():
            col = df[var].copy()
            col["estimate"] = col["estimate"] - theory[var]
            ax.errorbar(
                x="T",
                y="estimate",
                yerr="error",
                data=col.reset_index(drop=False),
                label=f"{sampler.value}",
                alpha=0.5,
            )
            ax.set_title(f"{var}")
    ax.legend()
    for i in range(2):
        axes[1, i].set_xlabel("T")
        axes[i, 0].set_ylabel(r"$\Delta$")
    return axes
