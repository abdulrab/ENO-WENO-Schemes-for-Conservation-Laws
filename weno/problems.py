"""Standard scalar test problems for the WENO examples."""

from __future__ import annotations

import numpy as np

from .core import Flux, Speed


def linear_advection(a: float = 1.0) -> tuple[Flux, Speed]:
    """Return flux and characteristic-speed functions for u_t + a u_x = 0."""

    def flux(u: np.ndarray) -> np.ndarray:
        return a * u

    def speed(u: np.ndarray) -> float:
        del u
        return abs(a)

    return flux, speed


def burgers() -> tuple[Flux, Speed]:
    """Return flux and characteristic-speed functions for inviscid Burgers."""

    def flux(u: np.ndarray) -> np.ndarray:
        return 0.5 * u**2

    def speed(u: np.ndarray) -> float:
        return float(np.max(np.abs(u)))

    return flux, speed
