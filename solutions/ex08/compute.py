import subprocess as proc
from typing import Tuple

import numpy as np
from numpy.typing import NDArray

from global_utils import executable, DEFAULT_PRIMES, DEFAULT_SEEDS
from ex08.paths import run_results, annealing_path


def run_simulated_annealing(run: int,
                            n_T_steps: int,
                            n_explore_steps: int,
                            T_bounds: Tuple[float, float],
                            n_blocks: int,
                            block_size: int,
                            initial_guess: Tuple[float, float],
                            sample_stddev: float,
                            use_cached: bool,
                            p_line: int = 0) -> proc.CompletedProcess[bytes]:
    """Searches for the pair (µ, σ) which minimizes the system's hamiltonian using the Simulated Annealing algorithm; temperature is reduced using uniform-logarithmic steps.
    Parameters are sampled by the Metropolis algorithm using a gaussian trasition with defined standard deviation.

    Args:
        n_T_steps (int): Number of temperature adjustments.
        n_explore_steps (int): Number of steps taken while keeping the temperature constant.
        T_bounds (Tuple[float, float]): Pair (initial temperature, final temperature).
        n_blocks (int): Number of blocks used to estimate <H> every time a pair of parameter is sampled.
        block_size (int): Block size.
        initial_guess (Tuple[float, float]): Initial parameters guess.
        sample_stddev (float): Standard deviation used by the normal sampler which proposes the next parameters candidates.
        p_line (int, optional): Line in the primes file.

    Returns:
        CompletedProcess[bytes]: Process state.
    """
    if use_cached and annealing_path(run).exists():
        return
    return proc.run([executable("08_2"),
                     f"--out={run_results(run)}",
                     f"--primes_path={DEFAULT_PRIMES}",
                     f"--seeds_path={DEFAULT_SEEDS}",
                     f"--primes_line={p_line}",
                     f"--n_steps={n_T_steps}",
                     f"--n_explore={n_explore_steps}",
                     f"--n_blocks={n_blocks}",
                     f"--block_size={block_size}",
                     f"--T0={T_bounds[0]}",
                     f"--Tf={T_bounds[1]}",
                     f"--p0={initial_guess[0]},{initial_guess[1]}",
                     f"--stddev={sample_stddev}",
                     ])


def run_estimation(run: int,
                   n_blocks: int,
                   block_size: int,
                   sampling_bound: float,
                   n_bins: int,
                   mu: float, sigma: float,
                   p_line: int = 0):
    """Provides an estimation of <H>_ψ_{µ,σ} and computes an histogram of |ψ_{µ,σ}|²

    Args:
        n_blocks (int): Number of blocks used for the estimation of <H>
        block_size (int): Block size. Must be odd.
        sampling_bound (float): x will be sampled in (-sampling_bound, +sampling_bound).
        n_bins (int): Number of bins in ψ's histogram.
        mu (float): ψ's µ.
        sigma (float): ψ's σ.
        p_line (int, optional): Line in the primes file. Defaults to 0.

    Returns:
        Any: Process state.
    """
    if n_bins % 2 != 1:
        print("The number of bins should be odd, so that the following computation can normalize on psi(0)")
    return proc.run([executable("08_2_psi"),
                     f"--out={run_results(run)}",
                     f"--primes_path={DEFAULT_PRIMES}",
                     f"--seeds_path={DEFAULT_SEEDS}",
                     f"--primes_line={p_line}",
                     f"--n_blocks={n_blocks}",
                     f"--block_size={block_size}",
                     f"--bounds={-sampling_bound},{sampling_bound}",
                     f"--n_bins={n_bins}",
                     f"--mu={mu}",
                     f"--sigma={sigma}"
                     ])


def compute_states(n_iter: int, sample_bounds: Tuple[float, float]) -> Tuple[NDArray, NDArray, NDArray]:
    """Estimates the energy levels using the matrix method.

    Returns:
        Tuple[NDArray, NDArray, NDArray]: Tuple (energies, state points, state fn)
    """
    def Vpot(x):
        return (x**2 - 2.5)*x**2

    hbar = 1
    m = 1

    # Step sizes
    x = np.linspace(*sample_bounds, n_iter)
    dx = x[1] - x[0]  # the step size
    V = Vpot(x)

    # The central differences method: f" = (f_1 - 2*f_0 + f_-1)/dx^2

    CDiff = np.diag(np.ones(n_iter-1), -1)-2 * \
        np.diag(np.ones(n_iter), 0)+np.diag(np.ones(n_iter-1), 1)
    # np.diag(np.array,k) construct a "diagonal" matrix using the np.array
    # The default is k=0. Use k>0 for diagonals above the main diagonal,
    # and k<0 for diagonals below the main diagonal

    # Hamiltonian matrix
    H = (-(hbar**2)*CDiff)/(2*m*dx**2) + np.diag(V)

    # Compute eigenvectors and their eigenvalues
    E, psi = np.linalg.eigh(H)

    # Take the transpose & normalize
    psi = np.transpose(psi)
    psi = psi/np.sqrt(dx)

    return E, x, psi
