from typing import Optional
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

ROOT = Path("..")
DEFAULT_PRIMES = ROOT / "data" / "primes" / "Primes"
DEFAULT_SEEDS = ROOT / "data" / "seeds" / "seed.in"


def results_dir(section: str) -> Path:
    """Path to the directory containing section's results

    Args:
        section (str): The section's name.
    """
    return ROOT / "results" / section


def executable(exercise: str) -> Path:
    """Path to the executable named "exercise"

    Args:
        exercise (str): The executable name.
    """
    return ROOT / "bin" / exercise


def default_output(exercise: str) -> Path:
    """A common scheme for the exercises' outputs

    Args:
        exercise (str): The exercise's number.
    """
    section = exercise.split("_")[0]
    return results_dir(section) / (exercise + ".csv")


def plot_trajectory(x, y, ax: plt.Axes, data: Optional[pd.DataFrame] = None, **kwargs):
    if data is not None:
        x = data[x]
        y = data[y]
    c = np.arange(len(x)) / len(x)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # norm = plt.Normalize(c[0], c[1])
    lc = LineCollection(segments, **kwargs)
    lc.set_array(c)
    ax.add_collection(lc)
