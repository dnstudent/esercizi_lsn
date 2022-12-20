from typing import Tuple, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# Autocorrelation function implemented using fft
from statsmodels.tsa.stattools import acf


def plot_autocorrs(data: Dict[str, pd.Series], nlags: int = 300, **kwargs) -> Tuple[plt.Figure, plt.Axes]:
    """Plots the autocorrelation functions of the given phases' measures.

    Args:
        phases (Dict[str, pd.Series]): Dictionary phase -> measures.
        nlags (int, optional): Number of lags for the autocorrelation fn.
        **kwargs: Keyworded arguments passed to plt.subplots.

    Returns:
        Tuple[plt.Figure, plt.Axes]: Matplotlib Figure and Axes
    """
    fig, ax = plt.subplots(1, 1, constrained_layout=True, **kwargs)
    acfns = []
    for phase, df in data.items():
        acfns.append(acf(df, nlags=nlags, fft=True))
        ax.plot(acfns[-1], label=phase)
    ax.axhline(np.exp(-1), label="$e^{-1}$", c="red")
    ax.legend()
    fig.suptitle("Autocorrelation functions")
    return acfns, fig, ax
