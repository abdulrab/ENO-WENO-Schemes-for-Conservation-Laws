"""Smooth periodic linear-advection demonstration.

Run from the repository root:

    python examples/linear_advection.py
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from weno import l1_error, linear_advection, linf_error, solve_periodic


def initial_condition(x: np.ndarray) -> np.ndarray:
    return np.sin(2.0 * np.pi * x)


def exact_solution(x: np.ndarray, t: float, a: float = 1.0) -> np.ndarray:
    return np.sin(2.0 * np.pi * ((x - a * t) % 1.0))


def main() -> None:
    a = 1.0
    t_final = 1.0
    flux, speed = linear_advection(a)

    sol = solve_periodic(
        initial_condition,
        flux,
        speed,
        n=200,
        t_final=t_final,
        cfl=0.4,
    )

    exact = exact_solution(sol.x, t_final, a)
    l1 = l1_error(sol.u, exact, sol.dx)
    linf = linf_error(sol.u, exact)

    print(f"steps = {sol.steps}")
    print(f"L1 error   = {l1:.8e}")
    print(f"Linf error = {linf:.8e}")

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(sol.x, exact, label="Exact", linewidth=2.0)
    ax.plot(sol.x, sol.u, "--", label="WENO5 + SSP-RK3", linewidth=1.5)
    ax.set_xlabel("x")
    ax.set_ylabel("u(x,t)")
    ax.set_title("Linear advection at t = 1")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()

    out = Path(__file__).resolve().parents[1] / "figures" / "linear_advection.svg"
    fig.savefig(out)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
