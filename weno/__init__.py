"""Small WENO5 toolkit for one-dimensional scalar conservation laws."""

from .core import (
    Solution,
    numerical_flux,
    solve_periodic,
    spatial_operator,
    ssprk3_step,
)
from .diagnostics import l1_error, l2_error, linf_error, observed_order
from .problems import burgers, linear_advection

__all__ = [
    "Solution",
    "burgers",
    "l1_error",
    "l2_error",
    "linear_advection",
    "linf_error",
    "numerical_flux",
    "observed_order",
    "solve_periodic",
    "spatial_operator",
    "ssprk3_step",
]
