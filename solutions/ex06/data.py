from typing import Tuple, List
import pathlib

import pandas as pd

from global_utils import results_dir
from .vars import *


def samplers(fn):
    """Python decorator. Automates the process of declaring functions which take multiple samplers and return results depending on the sampler.

    Args:
        fn (function): A function which takes an argument labeled 'sampler'
    """

    def wrapper(*args, sampler, **kwargs):
        if not hasattr(sampler, "__iter__"):
            sampler = [sampler]
        return {s: fn(*args, sampler=s, **kwargs) for s in sampler}

    return wrapper


def data_dir(T: float) -> pathlib.Path:
    """The directory containing measures at temperature T

    Args:
        T (float): Temperature.
        sampler (Sampler): MCMC algorithm used to sample.

    Returns:
        pathlib.Path: The directory containing to the measures at temperature T
    """
    return results_dir(SECTION) / str(T)


def spins_path(T: float, h: float, sampler: Sampler) -> pathlib.Path:
    """The path to the last sampled state from a previous run,

    Args:
        T (float): Temperature.
        h (float): External magnetic field.
        sampler (Sampler): MCMC algorithm used to sample.

    Returns:
        pathlib.Path: The path to the sampled state.
    """
    return data_dir(T) / f"{sampler.value}_{h:06f}_spins.csv"


def equilibration_path(T: float, h: float, n: int, sampler: Sampler) -> pathlib.Path:
    """Path to the equilibration results

    Args:
        T (float): Temperature.
        h (float): _description_
        n (int): _description_
        sampler (Sampler): _description_
    """
    return data_dir(T) / f"{sampler.value}_{h:06f}_warmup{n}.csv"


def autocorr_path(T: float, h: float, sampler: Sampler) -> pathlib.Path:
    """Path to the autocorrelation results

    Args:
        T (float): Temperature.
        h (float): _description_
        sampler (Sampler): _description_
    """
    return data_dir(T) / f"{sampler.value}_{h:06f}_autocorrelation.csv"


@samplers
def read_equilibrations(
    T: float, h: float, /, sampler: Sampler, nrows=None, skiprows=None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Reads the equilibration products

    Args:
        T (float): Temperature.
        h (float): External magnetic field, when active.
        sampler (Sampler): The MCMC algorithm used to sample points.

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Two separate runs of equilibration.
    """

    def eq_df(h, n, sampler: Sampler):
        return pd.read_csv(
            equilibration_path(T, h, n, sampler), nrows=nrows, skiprows=skiprows
        )

    dfs = []
    for n in [1, 2]:
        df0 = eq_df(0.0, n, sampler)
        df0["Sum_s"] = pd.NA
        df0.update(eq_df(h, n, sampler), overwrite=False)
        dfs.append(df0)
    return tuple(dfs)


@samplers
def read_autocorrelation(T: float, h: float, /, sampler: Sampler) -> pd.DataFrame:
    """Reads the autocorrelation performed on a previous warmup run.

    Args:
        T (float): _description_
        h (float): _description_
        sampler (Sampler): _description_
    """

    def autoc_df(h):
        return pd.read_csv(autocorr_path(T, h, sampler))

    df0 = autoc_df(0.0)
    df0["Sum_s"] = pd.NA
    df0.update(autoc_df(h), overwrite=False)
    return df0


@samplers
def read_state(T: float, h: float, /, sampler: Sampler) -> pd.DataFrame:
    """Reads the csv file containing the final state sampled by the MCMC algorithm for a system at temperature T during a previous run.

    Args:
        T (float): Temperature.
        h (float): External magnetic field.
        sampler (Sampler): MCMC algorithm used to sample.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the sampled state.
    """
    return pd.read_csv(spins_path(T, h, sampler))


@samplers
def read_measures(
    T: float, h: float, /, sampler: Sampler, nrows=None, skiprows=None
) -> pd.DataFrame:
    """Reads the csv file storing the thermodynamic variables estimations computed using the MCMC algorithm for a system at temperature T during a previous run.

    Args:
        T (float): Temperature.
        h (float): External magnetic field.
        sampler (Sampler): MCMC algorithm used to sample.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the variables estimations.
    """

    def file_path(h):
        return data_dir(T) / f"{sampler.value}_{h:06f}_variables.csv"

    df = pd.read_csv(file_path(0.0), nrows=nrows, skiprows=skiprows).join(
        pd.read_csv(file_path(h), nrows=nrows, skiprows=skiprows)
    )
    columns = pd.MultiIndex.from_tuples(
        map(tuple, df.columns.str.split("_", n=1)))
    df.columns = columns
    return df


def read_estimates(
    Ts: List[float], h: float, samplers: Union[Sampler, List[Sampler]]
) -> Dict[Sampler, pd.DataFrame]:
    if type(samplers) is Sampler:
        samplers = [samplers]
    last_estimates = {sampler: pd.DataFrame() for sampler in samplers}
    for T in Ts:
        estimations = read_measures(T, h, sampler=samplers)
        for sampler, df in estimations.items():
            l = df.tail(1)
            l.index = pd.Index([T], name="T")
            last_estimates[sampler] = pd.concat([last_estimates[sampler], l])
    return last_estimates


def read_results(
    Ts: List[float], h: float, samplers: Union[Sampler, List[Sampler]]
) -> Dict[Sampler, pd.DataFrame]:
    """Retrieves the estimation results and errors in a table for each sampler.

    Args:
        Ts (List[float]): Temperatures.
        h (float): External magnetic field.
        samplers (Union[Sampler, List[Sampler]]): MCMC samplers.

    Returns:
        Dict[Sampler, pd.DataFrame]: Results as a dictionary sampler->table.
    """
    return {
        sampler: df.rename_axis(index="T", columns=["Variable", "qty"])
        for sampler, df in read_estimates(Ts, h, samplers).items()
    }
