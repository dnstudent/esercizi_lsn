import pandas as pd

from .paths import results_path_from, distances_per_iter_path_from, TSP, Algo


"""Functions to return the results in a pandas DataFrame
"""


def data_from(tsp: TSP, algo: Algo, p_line: int) -> pd.DataFrame:
    df = pd.read_csv(results_path_from(tsp, algo, p_line))
    df.columns = list(range(len(df.columns)-1)) + [df.columns[-1]]
    return df


def distances_per_iter_from(tsp: TSP, algo: Algo, p_line: int) -> pd.DataFrame:
    return pd.read_csv(distances_per_iter_path_from(tsp, algo, p_line))
