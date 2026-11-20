# Navier-Stokes SIM

2D incompressible Navier-Stokes solver written using a staggered MAC grid and real-time visualization with raylib.

$$\frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u}\cdot\nabla)\mathbf{u} = -\frac{1}{\rho}\nabla p + \nu\nabla^2\mathbf{u} + \mathbf{f}, \qquad \nabla\cdot\mathbf{u}=0$$

## Overview

The goal of this project is to simulate a 2D channel flow past a circular cylinder and reproduce the Von Karman vortex street across different Reynolds numbers ($Re$). 

Each term of the equation (advection, viscosity, pressure projection, external forces) can be toggled independently at runtime to observe its effect on the flow.

## Numerical Scheme

The solver uses Chorin's projection method (operator splitting) on a Marker-and-Cell (MAC) staggered grid:

`add forces -> advect -> diffuse -> project`

* **Grid:** Staggered MAC layout (`CellField` at cell centers for $p$, $\text{div}$, dye; `FaceX` and `FaceY` at cell faces for $u$ and $v$) with ghost cells. Avoids checkerboard pressure modes and ensures exact discrete divergence-gradient consistency.
* **Advection:** Semi-Lagrangian backtracing with bilinear interpolation.
* **Diffusion:** Implicit viscosity solved via Gauss-Seidel relaxation.
* **Projection:** 5-point Poisson solver ($\nabla^2 p = \frac{1}{\Delta t}\nabla\cdot\mathbf{u}^*$) using Gauss-Seidel with warm start.

## Milestones

### Phase 1: Core Solver (Console)

- [x] **M1** — MAC staggered grid and data structures (`CellField`, `FaceX`, `FaceY`)
- [x] **M2** — Cross-grid interpolation and semi-Lagrangian advection
- [x] **M3** — Discrete divergence operator on MAC grid
- [x] **M4** — Poisson solver (Jacobi and Gauss-Seidel) with analytical verification
- [x] **M5** — Chorin pressure projection

### Phase 2: Visualization and Channel Physics

- [x] **M6** — Interactive window with `raylib` (dye and mouse forces)
- [x] **M7** — Implicit viscosity and live Reynolds number control
- [x] **M8** — Vorticity field visualization ($\omega = \nabla\times\mathbf{u}$)
- [x] **M9** — Rectangular channel boundary conditions and Poiseuille flow validation
- [x] **M10** — Cylinder obstacle masking and steady recirculation ($Re \approx 20$)
- [ ] **M11** — Von Karman vortex street and Strouhal number ($St$) measurement
- [ ] **M12** — Implement the capacity to activate/deactivate Navier-Stokes equation terms

## Dependencies

* C++17
* [raylib](https://www.raylib.com/) (from M6 onwards)
* Python 3 (`numpy`, `matplotlib`) for validation plots

## Build

* Execute build.sh
