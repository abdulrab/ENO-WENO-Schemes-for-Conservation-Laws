import numpy as np

from weno import (
    burgers,
    l1_error,
    linear_advection,
    observed_order,
    solve_periodic,
    spatial_operator,
)


def test_constant_state_is_preserved():
    flux, speed = burgers()
    u = np.full(64, 2.5)
    rhs = spatial_operator(u, dx=1.0 / u.size, flux=flux, max_characteristic_speed=speed)
    assert np.max(np.abs(rhs)) < 1.0e-12


def test_periodic_mass_is_conserved():
    flux, speed = burgers()

    def initial(x):
        return 0.5 + 0.25 * np.sin(2.0 * np.pi * x)

    n = 128
    dx = 1.0 / n
    x = dx * np.arange(n)
    mass0 = dx * np.sum(initial(x))

    sol = solve_periodic(initial, flux, speed, n=n, t_final=0.1, cfl=0.35)
    mass1 = sol.dx * np.sum(sol.u)

    assert abs(mass1 - mass0) < 5.0e-13


def test_linear_advection_refinement_reduces_error():
    flux, speed = linear_advection(1.0)

    def initial(x):
        return np.sin(2.0 * np.pi * x)

    errors = []
    for n in (40, 80):
        sol = solve_periodic(initial, flux, speed, n=n, t_final=0.25, cfl=0.4)
        exact = np.sin(2.0 * np.pi * ((sol.x - 0.25) % 1.0))
        errors.append(sol.dx * np.sum(np.abs(sol.u - exact)))

    assert errors[1] < 0.2 * errors[0]


def test_error_and_observed_order_utilities():
    numerical = np.array([1.0, 2.0, 3.0])
    exact = np.array([1.0, 1.0, 3.0])
    assert l1_error(numerical, exact, dx=0.5) == 0.5
    assert np.isclose(observed_order(8.0e-4, 1.0e-4), 3.0)
