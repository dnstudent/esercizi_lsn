from cmath import inf
import subprocess as proc
from typing import Any, List, Tuple, Optional

import numpy as np
from scipy.optimize import curve_fit, least_squares


from global_utils import executable, DEFAULT_PRIMES, DEFAULT_SEEDS
from .paths import zero_dir, measures_dir


def equilibrate_07_2(phase: str, p_line: int, use_cached: bool = False) -> Any:
    """Performs the equilibration as specified by the settings file ("input") in "zero_dir(phase)".

    Args:
        phase (str): Tag of the phase.
        p_line (int): Line in the prime numbers file.

    Returns:
        Any: Process state.
    """
    if not (use_cached and (measures_dir(2, phase, 'warmup')/'1'/'thermo.csv').exists()):
        proc.run([executable("07_2"),
                  f"--in={zero_dir(2, phase, 'warmup')}",
                  f"--out={measures_dir(2, phase, 'warmup')/'1'}",
                  f"--primes_path={DEFAULT_PRIMES}",
                  f"--seeds_path={DEFAULT_SEEDS}",
                  f"--primes_line={p_line}",
                  ], capture_output=True)

    if use_cached and (measures_dir(2, phase, 'warmup')/'2'/'thermo.csv').exists():
        return
    return proc.run([executable("07_2"),
                     f"--in={zero_dir(2, phase, 'warmup')}",
                     f"--out={measures_dir(2, phase, 'warmup')/'2'}",
                     f"--primes_path={DEFAULT_PRIMES}",
                     f"--seeds_path={DEFAULT_SEEDS}",
                     f"--primes_line={p_line+1}",
                     ], capture_output=True)


def run_07_2(phase: str, p_line: int, use_cached: bool = False) -> Any:
    """Performs the equilibration as specified by the settings file ("input") in "zero_dir(phase)".

    Args:
        phase (str): Tag of the phase.
        p_line (int): Line in the prime numbers file.

    Returns:
        Any: Process state.
    """
    if use_cached and (measures_dir(2, phase, 'run') / "thermo.csv").exists():
        return
    return proc.run([executable("07_2"),
                     f"--in={measures_dir(2, phase, 'warmup')/'1'}",
                     f"--out={measures_dir(2, phase, 'run')}",
                     f"--settings={zero_dir(2, phase, 'run') / 'input'}",
                     f"--primes_path={DEFAULT_PRIMES}",
                     f"--primes_line={p_line}",
                     ], capture_output=True)


def equilibrate_07_4(phase: str, p_line: int, samplers: List[str] = ["mc", "md"]) -> Any:
    """Performs the equilibration as specified by the settings file ("input") in "zero_dir(phase)".

    Args:
        phase (str): Tag of the phase.
        n_bins (int): Number of bins in the computation of the radial function.
        p_line (int): Line in the prime numbers file.
        samplers (List[str]): Samplers to use. Valid elements are "md" and "mc".
    Returns:
        Any: Process state.
    """
    samplers_flags = tuple([f"--{sampler}" for sampler in samplers])
    return proc.run([executable("07_4"),
                     f"--in_mc={zero_dir(4, phase, 'warmup', 'mc')}",
                     f"--in_md={zero_dir(4, phase, 'warmup', 'md')}",
                     f"--out={measures_dir(4, phase, 'warmup')}",
                     f"--primes_path={DEFAULT_PRIMES}",
                     f"--seeds_path={DEFAULT_SEEDS}",
                     f"--primes_line={p_line}",
                     f"--n_bins={0}",
                     "--warmup",
                     *samplers_flags
                     ], capture_output=True)


def run_07_4(phase: str, n_bins: int, p_line: int, samplers: List[str] = ["mc", "md"], start_from="warmup") -> Any:
    """Performs the equilibration as specified by the settings file ("input") in "zero_dir(phase)".

    Args:
        phase (str): Tag of the phase.
        n_bins (int): Number of bins in the computation of the radial function.
        p_line (int): Line in the prime numbers file.
        samplers (List[str]): Samplers to use. Valid elements are "md" and "mc".
    Returns:
        Any: Process state.
    """
    samplers_flags = tuple([f"--{sampler}" for sampler in samplers])
    return proc.run([executable("07_4"),
                     f"--in_mc={measures_dir(4, phase, start_from, 'mc')}",
                     f"--mc_settings={zero_dir(4, phase, 'run', 'mc') / 'input'}",
                     f"--in_md={measures_dir(4, phase, start_from, 'md')}",
                     f"--md_settings={zero_dir(4, phase, 'run', 'md') / 'input'}",
                     f"--out={measures_dir(4, phase, 'run')}",
                     f"--primes_path={DEFAULT_PRIMES}",
                     f"--seeds_path={DEFAULT_SEEDS}",
                     f"--primes_line={p_line}",
                     f"--n_bins={n_bins}",
                     *samplers_flags
                     ], capture_output=True)


def block_uncert(sample: np.ndarray, block_size: int) -> float:
    """Computes the block uncertainties of the quantity sampled in "sample".

    Args:
        sample (np.ndarray): Sample.
        block_size (int): Block size.

    Returns:
        float: The uncertainty estimation.
    """
    # trimming some values at the beginning so that the sample can be
    # partitioned in same-sized blocks.
    sample = np.array(sample[len(sample) % block_size:])
    # partitioning the sample through reshape, then the blocks' averages (A_i) are taken.
    block_avgs = sample.reshape((-1, block_size)).mean(axis=-1)
    # the averages' standard deviation / √n_blocks is the error.
    return block_avgs.std(ddof=1) / np.sqrt(block_avgs.shape[0])


def estimate_tc(acfn: np.ndarray,
                p0: Optional[float] = None,
                bounds: Tuple[float, float] = (0.0, np.inf),
                **kwargs) -> Tuple[float, float]:
    """Estimates the correlation time of a series given its autocorrelation function
    using the least squares method.

    Args:
        acfn (np.ndarray): Autocorrelation function.
        p0 (float): Initial guess.
        bounds (Tuple[float, float]): Estimation bounds.
        **kwargs: Arguments passed to scipy's curve_fit.

    Returns:
        float: Estimation of the autocorrelation time.
        float: Variance of the estimation.
    """
    def model(xs, t_c):
        return np.exp(-xs/t_c)

    if p0 is not None:
        p0 = np.array([p0])
    # curve_fit uses least_squares under the hood
    est, cov = curve_fit(model, np.arange(0, len(acfn)), acfn,
                         p0=p0, bounds=bounds, **kwargs)
    return est[0], cov[0, 0]
