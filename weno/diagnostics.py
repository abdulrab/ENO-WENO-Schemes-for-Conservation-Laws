"""Error norms and grid-refinement utilities."""

from __future__ import annotations

import numpy as np


def l1_error(numerical: np.ndarray, exact: np.ndarray, dx: float) -> float:
    """Discrete approximation to the L1 error integral."""
    numerical = np.asarray(numerical, dtype=float)
    exact = np.asarray(exact, dtype=float)
    return float(dx * np.sum(np.abs(numerical - exact)))


def l2_error(numerical: np.ndarray, exact: np.ndarray, dx: float) -> float:
    """Discrete approximation to the L2 error norm."""
    numerical = np.asarray(numerical, dtype=float)
    exact = np.asarray(exact, dtype=float)
    return float(np.sqrt(dx * np.sum((numerical - exact) ** 2)))


def linf_error(numerical: np.ndarray, exact: np.ndarray) -> float:
    """Discrete maximum-norm error."""
    numerical = np.asarray(numerical, dtype=float)
    exact = np.asarray(exact, dtype=float)
    return float(np.max(np.abs(numerical - exact)))


def observed_order(coarse_error: float, fine_error: float, refinement_ratio: float = 2.0) -> float:
    """Observed convergence order from two positive errors."""
    if coarse_error <= 0.0 or fine_error <= 0.0:
        raise ValueError("Errors must be positive.")
    if refinement_ratio <= 1.0:
        raise ValueError("refinement_ratio must be greater than 1.")
    return float(np.log(coarse_error / fine_error) / np.log(refinement_ratio))
