import subprocess as proc
from typing import Optional, Dict, Any
from pathlib import Path
import numpy as np

from global_utils import executable, DEFAULT_PRIMES, DEFAULT_SEEDS
from .paths import Algo, results_path_from, results_dir


def run_gp(input: Path, migration_length: int, n_migrations: int, pop_size: int, n_continents: int, mut_rate: float, algo: Algo, fusion_p: Optional[float], p_line: int = 0, use_cached: bool = True) -> Optional[proc.CompletedProcess[bytes]]:
    """Run the genetic process

    Args:
        input (Path): CSV of city coordinates.
        migration_length (int): Number of iterations between migrations.
        n_migrations (int): Number of migrations (when the MPI continents exchange data).
        pop_size (int): Total population size, which will be split between processes.
        n_processes (int): Number of MPI processes which will be launched.
        mut_rate (float): Mutation rate.
        algo (Algo): Crossover algorithm to use. Check the "Algo" enum.
        fusion_p (Optional[float]): "fusion" algorithm's probability to select "my2" against "exmod".
        p_line (int, optional): Line from which to take the prime number in "DEFAULT_PRIMES". Defaults to 0.
        use_cached (bool, optional): Whether to keep cached data, if found. Defaults to True.
    """
    if fusion_p is None or np.isnan(fusion_p):
        assert (algo is not Algo.FUSION)
        p = 0.0
    else:
        p = fusion_p
    if not (use_cached and results_path_from(algo, n_continents, fusion_p).exists()):
        return proc.run(["mpirun",
                         "-c", f"{n_continents}",
                         f"{executable('10_2')}",
                         f"--out={results_dir(n_continents, fusion_p)}",
                         f"--in={input}",
                         f"--crossover={algo.value}",
                         f"--migration_length={migration_length}",
                         f"--n_migrations={n_migrations}",
                         f"--pop_size={pop_size}",
                         f"--mut_rate={mut_rate}",
                         f"--fusion_p={p}",
                         f"--primes_path={DEFAULT_PRIMES}",
                         f"--seeds_path={DEFAULT_SEEDS}",
                         f"--primes_line={p_line}",
                         ])
