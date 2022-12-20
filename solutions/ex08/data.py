import pandas as pd

from .paths import min_energy_path, psi_path, annealing_path

# Read data from the default results paths as pandas DataFrame


def read_annealing(run: int) -> pd.DataFrame:
    return pd.read_csv(annealing_path(run))


def read_min_energy(run: int) -> pd.DataFrame:
    return pd.read_csv(min_energy_path(run))


def read_psi(run: int) -> pd.DataFrame:
    return pd.read_csv(psi_path(run))
