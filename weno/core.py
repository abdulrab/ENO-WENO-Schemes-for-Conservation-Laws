"""Fifth-order WENO finite-difference solver for 1-D scalar conservation laws.

The semidiscrete equation is

    u_t + f(u)_x = 0,

on a periodic uniform grid. Numerical fluxes use global Lax-Friedrichs
flux splitting, fifth-order Jiang-Shu WENO reconstruction, and SSP-RK3
time integration.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

Array = np.ndarray
Flux = Callable[[Array], Array]
Speed = Callable[[Array], float]


@dataclass(frozen=True)
class Solution:
    """Container returned by :func:`solve_periodic`."""

    x: Array
    u: Array
    t: float
    dx: float
    steps: int


def _weno5_left(v: Array, epsilon: float = 1.0e-6) -> Array:
    """Reconstruct v_{i+1/2} from the left on a periodic grid."""
    v = np.asarray(v, dtype=float)
    vm2 = np.roll(v, 2)
    vm1 = np.roll(v, 1)
    vp1 = np.roll(v, -1)
    vp2 = np.roll(v, -2)

    q0 = (2.0 * vm2 - 7.0 * vm1 + 11.0 * v) / 6.0
    q1 = (-vm1 + 5.0 * v + 2.0 * vp1) / 6.0
    q2 = (2.0 * v + 5.0 * vp1 - vp2) / 6.0

    beta0 = (
        (13.0 / 12.0) * (vm2 - 2.0 * vm1 + v) ** 2
        + 0.25 * (vm2 - 4.0 * vm1 + 3.0 * v) ** 2
    )
    beta1 = (
        (13.0 / 12.0) * (vm1 - 2.0 * v + vp1) ** 2
        + 0.25 * (vm1 - vp1) ** 2
    )
    beta2 = (
        (13.0 / 12.0) * (v - 2.0 * vp1 + vp2) ** 2
        + 0.25 * (3.0 * v - 4.0 * vp1 + vp2) ** 2
    )

    alpha0 = 0.1 / (epsilon + beta0) ** 2
    alpha1 = 0.6 / (epsilon + beta1) ** 2
    alpha2 = 0.3 / (epsilon + beta2) ** 2
    alpha_sum = alpha0 + alpha1 + alpha2

    return (alpha0 * q0 + alpha1 * q1 + alpha2 * q2) / alpha_sum


def _weno5_right(v: Array, epsilon: float = 1.0e-6) -> Array:
    """Reconstruct v_{i+1/2} from the right on a periodic grid."""
    v = np.asarray(v, dtype=float)
    vm1 = np.roll(v, 1)
    vp1 = np.roll(v, -1)
    vp2 = np.roll(v, -2)
    vp3 = np.roll(v, -3)

    q0 = (-vp3 + 5.0 * vp2 + 2.0 * vp1) / 6.0
    q1 = (2.0 * vp2 + 5.0 * vp1 - v) / 6.0
    q2 = (11.0 * vp1 - 7.0 * v + 2.0 * vm1) / 6.0

    beta0 = (
        (13.0 / 12.0) * (vp1 - 2.0 * vp2 + vp3) ** 2
        + 0.25 * (3.0 * vp1 - 4.0 * vp2 + vp3) ** 2
    )
    beta1 = (
        (13.0 / 12.0) * (v - 2.0 * vp1 + vp2) ** 2
        + 0.25 * (v - vp2) ** 2
    )
    beta2 = (
        (13.0 / 12.0) * (vm1 - 2.0 * v + vp1) ** 2
        + 0.25 * (vm1 - 4.0 * v + 3.0 * vp1) ** 2
    )

    alpha0 = 0.1 / (epsilon + beta0) ** 2
    alpha1 = 0.6 / (epsilon + beta1) ** 2
    alpha2 = 0.3 / (epsilon + beta2) ** 2
    alpha_sum = alpha0 + alpha1 + alpha2

    return (alpha0 * q0 + alpha1 * q1 + alpha2 * q2) / alpha_sum


def numerical_flux(
    u: Array,
    flux: Flux,
    max_characteristic_speed: Speed,
    epsilon: float = 1.0e-6,
) -> Array:
    """Return WENO5 numerical flux H_{i+1/2} at all interfaces."""
    f = np.asarray(flux(u), dtype=float)
    alpha = float(max_characteristic_speed(u))
    if alpha < 0.0:
        raise ValueError("Maximum characteristic speed must be nonnegative.")

    f_plus = 0.5 * (f + alpha * u)
    f_minus = 0.5 * (f - alpha * u)

    return _weno5_left(f_plus, epsilon) + _weno5_right(f_minus, epsilon)


def spatial_operator(
    u: Array,
    dx: float,
    flux: Flux,
    max_characteristic_speed: Speed,
    epsilon: float = 1.0e-6,
) -> Array:
    """Evaluate L(u) = -d f(u)/dx using conservative WENO5 fluxes."""
    interface_flux = numerical_flux(
        u,
        flux=flux,
        max_characteristic_speed=max_characteristic_speed,
        epsilon=epsilon,
    )
    return -(interface_flux - np.roll(interface_flux, 1)) / dx


def ssprk3_step(
    u: Array,
    dt: float,
    dx: float,
    flux: Flux,
    max_characteristic_speed: Speed,
    epsilon: float = 1.0e-6,
) -> Array:
    """Advance one time step with the third-order SSP Runge-Kutta method."""

    def L(v: Array) -> Array:
        return spatial_operator(
            v,
            dx=dx,
            flux=flux,
            max_characteristic_speed=max_characteristic_speed,
            epsilon=epsilon,
        )

    u1 = u + dt * L(u)
    u2 = 0.75 * u + 0.25 * (u1 + dt * L(u1))
    return (1.0 / 3.0) * u + (2.0 / 3.0) * (u2 + dt * L(u2))


def solve_periodic(
    initial_condition: Callable[[Array], Array],
    flux: Flux,
    max_characteristic_speed: Speed,
    *,
    x_left: float = 0.0,
    x_right: float = 1.0,
    n: int = 200,
    t_final: float = 1.0,
    cfl: float = 0.4,
    epsilon: float = 1.0e-6,
) -> Solution:
    """Solve a scalar conservation law on a periodic uniform grid.

    Grid points are x_j = x_left + j*dx for j=0,...,n-1.
    The time step is selected adaptively from

        dt = CFL * dx / max |f'(u)|.

    If the characteristic speed is numerically zero, the solution is
    returned unchanged at ``t_final``.
    """
    if n < 8:
        raise ValueError("WENO5 requires at least 8 periodic grid points.")
    if x_right <= x_left:
        raise ValueError("x_right must be greater than x_left.")
    if t_final < 0.0:
        raise ValueError("t_final must be nonnegative.")
    if not (0.0 < cfl <= 1.0):
        raise ValueError("cfl must satisfy 0 < cfl <= 1.")

    dx = (x_right - x_left) / n
    x = x_left + dx * np.arange(n, dtype=float)
    u = np.asarray(initial_condition(x), dtype=float).copy()
    if u.shape != x.shape:
        raise ValueError("initial_condition must return an array with shape (n,).")

    t = 0.0
    steps = 0
    time_tol = 100.0 * np.finfo(float).eps * max(1.0, t_final)

    while t < t_final - time_tol:
        speed = float(max_characteristic_speed(u))
        if speed < 0.0:
            raise ValueError("Maximum characteristic speed must be nonnegative.")
        if speed <= np.finfo(float).eps:
            t = t_final
            break

        dt = min(cfl * dx / speed, t_final - t)
        u = ssprk3_step(
            u,
            dt=dt,
            dx=dx,
            flux=flux,
            max_characteristic_speed=max_characteristic_speed,
            epsilon=epsilon,
        )
        t += dt
        steps += 1

    return Solution(x=x, u=u, t=t_final, dx=dx, steps=steps)
