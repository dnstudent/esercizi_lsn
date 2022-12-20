from enum import Enum
import pathlib
from typing import Optional

from numpy import isnan

from global_utils import results_dir


class Algo(Enum):
    EXERCISE = "ex"
    EXERCISEMOD = "exmod"
    MYALGO2 = "my2"
    FUSION = "fusion"
    DUMMY = "dummy"


DATA_PATH = results_dir("10")
CAPITALS_PATH = DATA_PATH / "capitals.csv"


def results_dir(n_continents: int, fusion_p: Optional[float]) -> pathlib.Path:
    if fusion_p is None or isnan(fusion_p):
        return DATA_PATH / f"{n_continents}"
    return DATA_PATH / f"{n_continents}" / f"fusion_{fusion_p:.3f}"


def results_path_from(algo: Algo, n_continents: int, fusion_p: Optional[float]) -> pathlib.Path:
    """Returns the path to the results of a TSP using a given algorithm.

    Args:
        algo (Algo): Which crossover algorithm was used.
        n_continents (int): Number of continents that was used to run the algorithm.
        fusion_p (Optional[float], optional): p parameter for the FUSION algorithm: it expresses the probability of 

    Returns:
        pathlib.Path: Path.
    """
    assert (algo is not Algo.FUSION or (
        algo is Algo.FUSION and not (fusion_p is None or isnan(fusion_p))))
    return results_dir(n_continents, fusion_p) / f"{algo.value}.csv"
