from enum import Enum
import pathlib

from global_utils import results_dir


class Algo(Enum):
    EXERCISE = "ex"
    EXERCISEMOD = "exmod"
    MYALGO2 = "my2"
    FUSION = "fusion"
    DUMMY = "dummy"


class TSP(Enum):
    CIRCLE = "circle"
    SQUARE = "square"
    AMERICA = "america"


DATA_PATH = results_dir("09")


def results_dir(tsp: TSP, p_line: int) -> pathlib.Path:
    return DATA_PATH / f"{tsp.value}" / f"{p_line}"


def results_path_from(tsp: TSP, algo: Algo, p_line: int) -> pathlib.Path:
    """Returns the path to the results of a TSP using a given algorithm.

    Args:
        tsp (TSP): Which TSP was run.
        algo (Algo): Which crossover algorithm was used.

    Returns:
        pathlib.Path: Path.
    """
    return results_dir(tsp, p_line) / f"{algo.value}.csv"


def distances_per_iter_path_from(tsp: TSP, algo: Algo, p_line: int) -> pathlib.Path:
    """Returns the path to the average distance per iteration for the better half of the population.

    Args:
        tsp (TSP): Which TSP was run.
        algo (Algo): Which crossover algorithm was used.

    Returns:
        pathlib.Path: Path.
    """
    return results_dir(tsp, p_line) / f"{algo.value}_stats.csv"
