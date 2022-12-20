from typing import Tuple, Callable

import numpy as np
from numpy.typing import ArrayLike


def measure(mean: ArrayLike, sigma: ArrayLike, rng: np.random.Generator) -> np.ndarray:
    """Simulates the measurement of a variable affected by gaussian error.

    Args:
        mean (ArrayLike): Measurement's true value.
        sigma (ArrayLike): Gaussian error.
        rng (np.random.Generator): Random number generator.

    Returns:
        np.ndarray: The measure, affected by error.
    """
    return rng.normal(mean, sigma)


def generate_data(domain: Tuple[ArrayLike, ArrayLike], f: Callable[[ArrayLike], ArrayLike], sigma: float, N: int, rng: np.random.Generator) -> Tuple[np.ndarray, np.ndarray]:
    """Uniformly generates a number of values in a domain and performs measures of a property affected by gaussian error on them.

    Args:
        domain (Tuple[ArrayLike, ArrayLike]): Independent variable's domain.
        f (Callable[[ArrayLike], ArrayLike]): Property to be measured.
        sigma (float): Gaussian error on the property's measurement.
        N (int): Number of values.
        rng (np.random.Generator): Numpy random number generator.

    Returns:
        Tuple[np.ndarray, np.ndarray]: Tuple (values, measures)
    """
    if len(np.shape(domain)) < 2:
        domain = np.expand_dims(domain, axis=-1)
    x = rng.uniform(domain[0], domain[1], (N,) + np.shape(domain)[1:])
    y = measure(f(x), sigma, rng)
    return x, y
