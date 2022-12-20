import subprocess as proc
from typing import Optional
from pathlib import Path
import os

import numpy as np
import pandas as pd

from global_utils import executable, DEFAULT_PRIMES, DEFAULT_SEEDS
from ex09.paths import Algo, TSP, results_path_from, results_dir


def generate_circle_coordinates(path: Path, radius: float = 1.0):
    ths = np.linspace(0, 2*np.pi, 34, endpoint=False)
    points = radius * np.array([np.cos(ths), np.sin(ths)])
    if not path.exists():
        os.makedirs(path.parent, exist_ok=True)
    pd.DataFrame(points).T.to_csv(path, index=False, header=False)


def generate_square_coordinates(path: Path, seed: int, side: float = 2.0):
    rng = np.random.default_rng(seed)
    points = rng.uniform(low=np.full(2, -side/2),
                         high=np.full(2, side/2), size=(34, 2))
    if not path.exists():
        os.makedirs(path.parent, exist_ok=True)
    pd.DataFrame(points).to_csv(path, index=False, header=False)


def run_ga(input: Path, n_iter: int, pop_size: int, mut_rate: float, algo: Algo, tsp: TSP, fusion_p: float = 0.4, p_line: int = 0, use_cached: bool = True) -> Optional[proc.CompletedProcess[bytes]]:
    """Run the genetic process

    Args:
        input (Path): CSV of city coordinates.
        n_iter (int): Number of iterations.
        pop_size (int): Population size.
        mut_rate (float): Mutation rate.
        algo (Algo): Crossover algorithm to use. Check the "Algo" enum.
        fusion_p (float, optional): "fusion" algorithm's probability to select "my2". Defaults to 0.4.
        p_line (int, optional): Line from which to take the prime number in "DEFAULT_PRIMES". Defaults to 0.
        use_cached (bool, optional): Whether to keep cached data, if found. Defaults to True.
    """
    if not (use_cached and results_path_from(tsp, algo, p_line).exists()):
        return proc.run([executable("09_1"),
                         f"--out={results_dir(tsp, p_line)}",
                         f"--in={input}",
                         f"--crossover={algo.value}",
                         f"--n_iter={n_iter}",
                         f"--pop_size={pop_size}",
                         f"--mut_rate={mut_rate}",
                         f"--fusion_p={fusion_p}",
                         f"--primes_path={DEFAULT_PRIMES}",
                         f"--seeds_path={DEFAULT_SEEDS}",
                         f"--primes_line={p_line}",
                         ])
