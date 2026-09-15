# ENO/WENO Schemes for Conservation Laws

This repository contains an academic presentation on **Essentially Non-Oscillatory (ENO)** and **Weighted Essentially Non-Oscillatory (WENO)** methods for hyperbolic conservation laws and Riemann problems.

The project was prepared at **Sukkur IBA University** as part of my graduate work in numerical analysis and computational methods for partial differential equations.

## Project status

This repository currently contains the **theoretical/proposal-stage presentation** rather than a complete numerical implementation. The existing material introduces the finite-volume framework, ENO reconstruction, stencil selection, numerical fluxes, and the method-of-lines viewpoint.

A full computational implementation of ENO/WENO methods is a natural next stage of the project.

## Mathematical setting

The presentation considers one-dimensional hyperbolic conservation laws of the form

\[
\frac{\partial u}{\partial t} + \frac{\partial f(u)}{\partial x} = 0,
\]

where \(u(x,t)\) is the conserved quantity and \(f(u)\) is the corresponding flux.

For a finite-volume discretization, the evolution of cell averages is written schematically as

\[
\frac{d\bar u_j}{dt}
+ \frac{1}{\Delta x}
\left(
\hat f_{j+1/2}-\hat f_{j-1/2}
\right)=0.
\]

The main numerical challenge is to achieve high-order accuracy in smooth regions while avoiding spurious oscillations near discontinuities such as shocks.

## Topics covered

- Hyperbolic partial differential equations
- Conservation laws
- Riemann problems
- Finite-volume methods
- ENO reconstruction
- Adaptive stencil selection using divided differences
- Interface reconstruction
- Numerical fluxes
- Method of lines
- Runge–Kutta time integration
- Motivation for WENO methods

## ENO reconstruction idea

ENO methods construct high-order polynomial reconstructions while adaptively choosing stencils that avoid regions with large solution variation.

A typical workflow is:

1. Start from cell averages \(\bar u_j\).
2. Select a stencil using divided differences.
3. Construct a polynomial whose cell averages agree with the selected neighboring cells.
4. Evaluate the polynomial at cell interfaces to obtain left and right reconstructed states.
5. Use these states in a numerical flux.
6. Advance the semi-discrete system in time, typically with a Runge–Kutta method.

This adaptive stencil selection is designed to reduce the oscillations that standard high-order polynomial approximations can generate near discontinuities.

## Repository contents

```text
.
├── main.tex          # Beamer presentation source
├── references.bib    # Bibliographic references
└── README.md         # Project overview
```

## References currently used

The presentation includes references to foundational and review literature, including work by:

- Chi-Wang Shu on ENO and WENO schemes
- Harten, Engquist, Osher, and Chakravarthy on high-order ENO methods
- Murillo and García-Navarro on Riemann problems for hyperbolic balance laws

See `references.bib` for the complete entries.

## Planned development

Future improvements may include:

- A reproducible Python implementation of ENO reconstruction
- WENO3 and WENO5 reconstruction
- Standard scalar test problems such as linear advection and Burgers' equation
- Riemann-problem test cases
- Comparison of numerical fluxes
- CFL and stability experiments
- Error and convergence studies
- Shock/rarefaction visualizations
- Extension toward nonlinear systems of conservation laws

## Related research interests

My broader research interests include numerical analysis, numerical PDEs, hyperbolic conservation laws, Riemann problems, finite-volume methods, two-phase flow, and scientific computing.

For more information, see my academic portfolio:

https://abdulrab.github.io/

## Author

**Abdul Rab**  
Lecturer in Mathematics  
Sukkur, Pakistan
