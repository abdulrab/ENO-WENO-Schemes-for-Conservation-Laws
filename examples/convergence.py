"""Grid-refinement study for smooth periodic linear advection.

Because the spatial discretization is WENO5 but time integration is SSP-RK3
with dt proportional to dx, sufficiently fine-grid results are expected to
approach third-order *fully discrete* convergence.

Run from the repository root:

    python examples/convergence.py
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from weno import l1_error, linear_advection, observed_order, solve_periodic


def initial_condition(x: np.ndarray) -> np.ndarray:
    return np.sin(2.0 * np.pi * x)


def exact_solution(x: np.ndarray, t: float) -> np.ndarray:
    return np.sin(2.0 * np.pi * ((x - t) % 1.0))


def main() -> None:
    flux, speed = linear_advection(1.0)
    ns = np.array([40, 80, 160, 320], dtype=int)
    errors = []

    for n in ns:
        sol = solve_periodic(
            initial_condition,
            flux,
            speed,
            n=int(n),
            t_final=1.0,
            cfl=0.4,
        )
        exact = exact_solution(sol.x, 1.0)
        error = l1_error(sol.u, exact, sol.dx)
        errors.append(error)

    errors = np.asarray(errors)
    orders = np.full_like(errors, np.nan)
    orders[1:] = [
        observed_order(errors[i - 1], errors[i])
        for i in range(1, len(errors))
    ]

    print("| N | L1 error | observed order |")
    print("|---:|---:|---:|")
    for n, err, order in zip(ns, errors, orders):
        order_text = "-" if np.isnan(order) else f"{order:.3f}"
        print(f"| {n} | {err:.8e} | {order_text} |")

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.loglog(ns, errors, "o-", label="Measured L1 error")
    reference = errors[-1] * (ns[-1] / ns) ** 3
    ax.loglog(ns, reference, "--", label="Third-order reference")
    ax.set_xlabel("Number of grid points N")
    ax.set_ylabel("L1 error")
    ax.set_title("Linear-advection grid refinement")
    ax.grid(alpha=0.25, which="both")
    ax.legend()
    fig.tight_layout()

    out = Path(__file__).resolve().parents[1] / "figures" / "convergence.svg"
    fig.savefig(out)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
