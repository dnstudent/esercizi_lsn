from enum import Enum
from typing import Dict, Union

import numpy as np

# Defining exercise-specific variables
SECTION = "06"
EQUILIBRATE = SECTION + "_warmup"
AUTOCORR = SECTION + "_correlation"
MEASURE = SECTION + "_measure"
PRIME_LINE = 1

VARIABLES = ["u", "c", "X", "m"]


class Sampler(Enum):
    Metropolis = "metropolis"
    Gibbs = "gibbs"


# Theoric values
def compute_U(
    J: float, beta: Union[float, np.ndarray], N: int
) -> Union[float, np.ndarray]:
    """Computes the internal energy per spin of a 1D, 50-spins, h=0 ising model"""
    th = np.tanh(J * beta)
    thN = th**N
    ch = 1 / th
    return -J * (th + ch * thN) / (1 + thN)


def compute_C(
    J: float, beta: Union[float, np.ndarray], N: int
) -> Union[float, np.ndarray]:
    """Computes the heat capacity per spin of a 1D, 50-spins, h=0 ising model"""
    th = np.tanh(J * beta)
    thN = th**N
    ch = 1 / th
    return ((beta * J) ** 2) * (
        ((1 + thN + (N - 1) * (th**2) + (N - 1) * (ch**2) * thN) / (1 + thN))
        - N * ((th + ch * thN) / (1 + thN)) ** 2
    )


def compute_X(
    J: float, beta: Union[float, np.ndarray], N: int
) -> Union[float, np.ndarray]:
    """Computes the magnetic susceptivity of a 1D, 50-spins, h=0 ising model"""
    th = np.tanh(J * beta)
    thN = th**N
    return beta * np.exp(2 * beta * J) * (1 - thN) / (1 + thN)


def compute_M(
    J: float, beta: Union[float, np.ndarray], h: float, N: int
) -> Union[float, np.ndarray]:
    """Computes the magnetization per spin of a 1D, 50-spins ising model with external magnetic field h"""
    bj = beta * J
    ebj = np.exp(bj)
    bh = beta * h
    l11 = ebj * np.cosh(bh)
    l12 = np.sqrt(np.square(ebj * np.cosh(bh)) - 2 * np.sinh(2 * bj))
    l1 = l11 + l12
    l2 = l11 - l12
    Z = np.power(l1, N) + np.power(l2, N)
    a = ebj * np.cosh(bh) / l12
    return (
        (np.power(l1, N - 1) * (1 + a) + np.power(l2, N - 1) * (1 - a))
        * ebj
        * np.sinh(bh)
        / Z
    )


def compute_theory(
    J: float, T: Union[float, np.ndarray], h: float, N: int
) -> Dict[str, Union[float, np.ndarray]]:
    """Computes the average values of U,C,X and M per spin of a 1D 50 spins ising model and stores them in a dictionary

    Args:
        J (float): Spin coupling.
        T (Union[float, np.ndarray]): Temperature(s).
        h (float): External magnetic field. Applied only when computing m.
        N (int): Number of spins.

    Returns:
        Dict[str, Union[float, np.ndarray]]: Dictionary {"variable": "theoric result(s)"}.
    """
    beta = 1 / T
    return {
        "u": compute_U(J, beta, N),
        "c": compute_C(J, beta, N),
        "X": compute_X(J, beta, N),
        "m": compute_M(J, beta, h, N),
    }
