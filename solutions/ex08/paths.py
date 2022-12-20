import pathlib
from global_utils import results_dir


def run_results(run: int) -> pathlib.Path:
    """Path to the results directory.

    Args:
        run (int): Run number
    """
    return results_dir("08") / str(run)


def annealing_path(run: int) -> pathlib.Path:
    """Path to the annealing results.

    Args:
        run (int): Run number
    """
    return run_results(run) / "annealing.csv"


def min_energy_path(run: int) -> pathlib.Path:
    """Path to the ground energy estimation.

    Args:
        run (int): Run number
    """
    return run_results(run) / "H_min.csv"


def psi_path(run: int) -> pathlib.Path:
    """Path to the approximated ground state sampling.

    Args:
        run (int): Run number
    """
    return run_results(run) / "psi.csv"
