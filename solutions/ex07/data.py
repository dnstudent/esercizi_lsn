import pandas as pd

from .paths import measures_dir


def read_equilibration_07_2(phase: str) -> pd.DataFrame:
    """Reads the equilibration data gathered for "phase".

    Args:
        phase (str): The phase tag.

    Returns:
        pd.DataFrame: A pandas DataFrame with the samples.
    """
    return (pd.read_csv(measures_dir(2, phase, 'warmup') / '1' / "thermo.csv"),
            pd.read_csv(measures_dir(2, phase, 'warmup') / '2' / "thermo.csv"))


def read_measures_07_2(phase: str) -> pd.DataFrame:
    """Reads the equilibration data gathered for "phase".

    Args:
        phase (str): The phase tag.

    Returns:
        pd.DataFrame: A pandas DataFrame with the samples.
    """
    return pd.read_csv(measures_dir(2, phase, 'run') / "thermo.csv")


def read_measures_07_4(phase: str, method: str) -> pd.DataFrame:
    """Reads the run data gathered for "phase".

    Args:
        phase (str): The phase tag.
        method (str): The method used (mc/md).

    Returns:
        pd.DataFrame: A pandas DataFrame with the samples.
    """
    try:
        return (pd.read_csv(measures_dir(4, phase, 'run', method) / "thermo.csv"),
                pd.read_csv(measures_dir(4, phase, 'run', method) / "g_r.csv"))
    except:
        return (pd.read_csv(measures_dir(4, phase, 'run', method) / "thermo.csv"), None)
