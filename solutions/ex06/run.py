from typing import List, Callable
import subprocess as proc
from unittest import skip

from statsmodels.tsa.stattools import acf

from global_utils import executable, DEFAULT_PRIMES, DEFAULT_SEEDS
from .data import *


def run_equilibration(
    J: float,
    h: float,
    n_spins: int,
    samplers: Union[Sampler, List[Sampler]],
    n_register: int,
    n_steps: int,
    save_spins: bool = True,
    use_cached: bool = False,
    prime_line: int = 1,
) -> Callable[[float], None]:
    """Returns a function in T which performs an equilibration of two Ising systems for the given number of steps
    registering the values of energy (H), sum of spins (Sum_s) and sum of product
    (Sum_s2).

    Args:
        J (float): Spin coupling.
        h (float): External magnetic field.
        n_spins (int): Number of spins.
        n_register (int): Number of steps for which variable values will be stored.
        samplers (Union[Sampler, List[Sampler]]): The MCMC sampler(s) to use.
        save_spins (bool, optional): Whether to save the final spin configuration. Defaults to True.

    Returns:
        Callable[[float], None]: A function of T.
    """
    assert n_steps >= n_register
    if type(samplers) is Sampler:
        samplers = [samplers]
    samplers_flag = [f"--{sampler.value}" for sampler in samplers]

    def run(T: float):
        if use_cached and len(list(data_dir(T).glob("*.csv"))) > 0:
            return

        return proc.run(
            [
                executable(EQUILIBRATE),
                f"--primes_path={DEFAULT_PRIMES}",
                f"--seeds_path={DEFAULT_SEEDS}",
                f"--primes_line={prime_line}",
                f"--n_steps={int(n_register)}",
                f"--block_size={1}",
                f"--n_warmup={int(n_steps - n_register)}",
                f"--out={data_dir(T)}",
                *samplers_flag,
                "--save_spins" if save_spins else "",
                f"--n_spins={n_spins}",
                f"--coupling={J}",
                f"--external_field={h}",
                f"--temperature={T}",
            ],
        )

    return run


def run_autocorr(
    h: float,
    samplers: Union[Sampler, List[Sampler]],
    n_lags: int,
    skip_steps: int,
    fft=True,
    use_cached=False
) -> Callable[[float], None]:
    """Returns a function in T which performs an autocorrelation analysis on a subset of the last run results.

    Args:
        T (float): Temperature.
        h (float): External magnetic field (when present).
        n_lags (int): Subset size.
        skip_steps (int): Number of steps at the beginning of the data which will not be included in the analysis.
        samplers (Union[Sampler, List[Sampler]]): Which sampler(s) to use.

    Returns:
        Callable[[float], None]: A function of T.
    """

    if type(samplers) is Sampler:
        samplers = [samplers]
    hs = [0.0, h]

    if fft:
        def run(T: float):
            eq_dfs = read_equilibrations(T, h, sampler=samplers)
            for sampler, (df, _) in eq_dfs.items():
                cfs_0 = df.iloc[skip_steps:].apply(acf, nlags=n_lags, fft=True)
                cfs_h = cfs_0.copy()
                cfs_0["Sum_s"] = 0.0
                cfs_0.to_csv(autocorr_path(T, 0, sampler), index=False)
                cfs_h["H"] = 0.0
                cfs_h["Sum_s2"] = 0.0
                cfs_h.to_csv(autocorr_path(T, h, sampler), index=False)

        return run

    def run(T: float):
        data_paths = {
            sampler: {h: equilibration_path(T, h, 1, sampler) for h in hs}
            for sampler in samplers
        }

        for sampler, path in data_paths.items():
            for actual_h in hs:
                if not (use_cached and autocorr_path(T, actual_h, sampler).exists()):
                    proc.run(
                        [
                            executable(AUTOCORR),
                            f"--in={path[actual_h]}",
                            f"--out={autocorr_path(T, actual_h, sampler)}",
                            f"--n_lags={n_lags}",
                            f"--skip={skip_steps}",
                        ]
                    )
    return run


def run_measures(
    J: float,
    h: float,
    n_spins: int,
    samplers: Union[Sampler, List[Sampler]],
    n_steps: int,
    block_size: int,
    save_spins: bool = False,
    prime_line: int = 1,
    use_cached=False,
) -> Callable[[float], None]:
    """Builds a function of T which performs a run of the experiment.

    Args:
        J (float): Spin coupling.
        h (float): External magnetic field. Applied only when measuring m.
        n_steps (int): Number of MC steps to take.
        block_size (int): The block size to use with the blocking method.
        samplers (Union[Sampler, List[Sampler]]): The MCMC sampler(s) to use.
        save_spins (bool, optional): Whether to save the final spin configuration. Defaults to True.
        resume (bool, optional): Whether to resume from the previous run. Defaults to False.
        prime_line (int, optional): Which line to use from the primes file. Defaults to 1.

    Returns:
        Callable[[float], None]: A function of T.
    """
    if not hasattr(samplers, "__iter__"):
        samplers = [samplers]
    samplers_flag = [f"--{sampler.value}" for sampler in samplers]

    def run(T: float):
        if not spins_path(T, h, samplers[0]).exists():
            resume = False
        else:
            resume = True
        if use_cached and not resume:
            return
        return proc.run(
            [
                executable(MEASURE),
                f"--primes_path={DEFAULT_PRIMES}",
                f"--seeds_path={DEFAULT_SEEDS}",
                f"--primes_line={prime_line}",
                f"--n_steps={n_steps}",
                f"--block_size={block_size}",
                f"--n_warmup=0",
                f"--out={data_dir(T)}",
                *samplers_flag,
                "--save_spins" if save_spins else "",
                "--resume" if resume else "",
                f"--n_spins={n_spins}",
                f"--coupling={J}",
                f"--external_field={h}",
                f"--temperature={T}",
            ]
        )

    return run


def run_measures2(
    J: float,
    h: float,
    n_spins: int,
    n_blocks: int,
    save_spins: bool = False,
    prime_line: int = 1,
) -> Callable[[float], None]:
    """Builds a function of T which performs a run of the experiment.

    Args:
        J (float): Spin coupling.
        h (float): External magnetic field. Applied only when measuring m.
        n_steps (int): Number of MC steps to take.
        block_size (int): The block size to use with the blocking method.
        samplers (Union[Sampler, List[Sampler]]): The MCMC sampler(s) to use.
        save_spins (bool, optional): Whether to save the final spin configuration. Defaults to True.
        resume (bool, optional): Whether to resume from the previous run. Defaults to False.
        prime_line (int, optional): Which line to use from the primes file. Defaults to 1.

    Returns:
        Callable[[float], None]: A function of T.
    """

    def run(args: Tuple[float, Sampler, int]):
        T, sampler, block_size = args
        if not spins_path(T, h, sampler).exists():
            resume = False
        else:
            resume = True
        return proc.run(
            [
                executable(MEASURE),
                f"--primes_path={DEFAULT_PRIMES}",
                f"--seeds_path={DEFAULT_SEEDS}",
                f"--primes_line={prime_line}",
                f"--n_steps={n_blocks*block_size}",
                f"--block_size={block_size}",
                f"--n_warmup=0",
                f"--out={data_dir(T)}",
                f"--{sampler.value}",
                "--save_spins" if save_spins else "",
                "--resume" if resume else "",
                f"--n_spins={n_spins}",
                f"--coupling={J}",
                f"--external_field={h}",
                f"--temperature={T}",
            ]
        )

    return run
