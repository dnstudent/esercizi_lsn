from typing import Optional, Dict, Any, List
import pandas as pd

from .paths import results_path_from, Algo


"""Functions to return the results in a pandas DataFrame
"""


def data_from(algo: Algo, n_continents: int, fusion_p: Optional[float]) -> pd.DataFrame:
    df = pd.read_csv(results_path_from(algo, n_continents, fusion_p))
    df.columns = list(range(len(df.columns)-1)) + [df.columns[-1]]
    return df


def collect_best_result_from(algo: Algo, n_continents: int, fusion_p: Optional[float], **kwargs):
    results = data_from(algo, n_continents,
                        fusion_p).sort_values("total_distance")
    best_individual = results.drop(columns="total_distance").iloc[0]
    shortest_distance = results["total_distance"].iloc[0]
    return {"shortest_distance": shortest_distance, "best_individual": best_individual}


def collect_best_results(configs: pd.DataFrame):
    return configs.join(configs.apply(lambda row: collect_best_result_from(**row), axis=1, result_type="expand")).sort_values("shortest_distance")
