"""Inviscid Burgers demonstration showing shock formation from smooth data.

Run from the repository root:

    python examples/burgers.py
"""

from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from weno import burgers, solve_periodic


def initial_condition(x: np.ndarray) -> np.ndarray:
    return 0.5 + 0.5 * np.sin(2.0 * np.pi * x)


def main() -> None:
    flux, speed = burgers()
    t_final = 0.5

    sol = solve_periodic(
        initial_condition,
        flux,
        speed,
        n=400,
        t_final=t_final,
        cfl=0.35,
    )

    initial = initial_condition(sol.x)
    mass_initial = sol.dx * np.sum(initial)
    mass_final = sol.dx * np.sum(sol.u)

    print(f"steps = {sol.steps}")
    print(f"initial mass = {mass_initial:.12e}")
    print(f"final mass   = {mass_final:.12e}")
    print(f"mass change  = {mass_final - mass_initial:.3e}")

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(sol.x, initial, label="Initial condition", linewidth=1.5)
    ax.plot(sol.x, sol.u, label=f"WENO5 solution, t={t_final:g}", linewidth=2.0)
    ax.set_xlabel("x")
    ax.set_ylabel("u(x,t)")
    ax.set_title("Inviscid Burgers equation: nonlinear steepening")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()

    out = Path(__file__).resolve().parents[1] / "figures" / "burgers.svg"
    fig.savefig(out)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
