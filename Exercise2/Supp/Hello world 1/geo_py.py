import numpy as np


def geo(z: float, n: int) -> float:
    """A geometric series up to n. sum_up_to_n(z^n)"""
    return sum([z**i for i in range(n)])


def geo_np(z: float, n: int):
    """A geometric series up to n. sum_up_to_n(z^n)"""
    return np.sum([np.full((n,), z) ** np.arange(n)])
