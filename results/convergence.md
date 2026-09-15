# Numerical verification

These results were generated from the code in this repository using the committed example scripts.

## Linear advection: grid refinement

Problem:

\[
u_t + u_x = 0,
\qquad
u(x,0)=\sin(2\pi x),
\qquad x\in[0,1),
\]

with periodic boundary conditions and final time \(t=1\).

Parameters:

- WENO5 spatial reconstruction;
- global Lax-Friedrichs flux splitting;
- SSP-RK3 time stepping;
- CFL = 0.4.

| N | L1 error | observed order |
|---:|---:|---:|
| 40 | 8.58907845e-05 | - |
| 80 | 6.56116358e-06 | 3.710 |
| 160 | 6.89533105e-07 | 3.250 |
| 320 | 8.21050906e-08 | 3.070 |

The spatial reconstruction is formally fifth order in smooth regions. Because the SSP-RK3 time step is chosen proportional to \(\Delta x\), the fully discrete asymptotic order is limited by the third-order temporal discretization. The observed order approaches approximately three on finer meshes.

## Linear advection: N = 200

At \(t=1\):

- L1 error = `3.45007299e-07`
- Linf error = `5.43318761e-07`
- time steps = `500`

## Inviscid Burgers equation

Problem:

\[
u_t+\left(\frac{u^2}{2}\right)_x=0,
\qquad
u(x,0)=\frac12+\frac12\sin(2\pi x).
\]

Parameters:

- N = 400;
- final time = 0.5;
- CFL = 0.35.

Measured mass:

- initial mass = `5.000000000000e-01`
- final mass = `5.000000000000e-01`
- change = `-1.721e-14`

The near-machine-precision difference is consistent with the conservative periodic flux-difference formulation.

## Reproduce

Run:

```bash
python examples/linear_advection.py
python examples/burgers.py
python examples/convergence.py
```

and verify the implementation with:

```bash
python -m pytest -q
```
