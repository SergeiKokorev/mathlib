from typing import Tuple


"""
This module provides some helpful mathematical functions
"""

def factorial(n: int) -> int:
    if n == 0:
        return 1
    elif n < 0:
        return -1
    else:
        return n * factorial(n - 1)


def linspace(start: float, stop: float, num_points: int) -> list:
    delta = (stop - start) / (num_points - 1)
    return [start + i * delta for i in range(num_points)]


def arange(start: float=0.0, stop: float=1.0, step: float=0.1):
    np = round((stop - start) / step)
    return [start + i * step for i in range(np)]
