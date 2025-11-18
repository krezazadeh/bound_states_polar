# Polar Gravitational Bound States Inside a Schwarzschild Black Hole

This program computes gravitational bound-state eigenfrequencies for polar perturbations inside a Schwarzschild black hole.

The eigenfrequencies are purely imaginary (omega = i*omegaI) and correspond to regular, normalizable solutions of the polar master equation.

The program uses a random-walk minimization method, Wronskian matching, and MPI parallelization to determine the value of omegaI that satisfies the internal boundary conditions.

## Physical Background

Polar gravitational perturbations inside a Schwarzschild black hole satisfy a second-order ODE with singular points at:
  - r = 0 (center)
  - r = 1 (horizon, with 2GM = 1)

To construct regular solutions:
  1. A power-series expansion (Frobenius method) is used near r = 0.
  2. The ODE is integrated outward (“minus” solution).
  3. Another solution is integrated inward from the horizon (“plus” solution).
  4. The Wronskian W = Zminus * Zplus' − Zplus * Zminus' is evaluated.
The correct eigenvalue is the one that minimizes |W|.

## Numerical Method

- The ODE system is solved using an eighth-order Runge–Kutta (RK8) integrator.
- A power-series expansion of order n_series provides initial conditions.
- A random-search minimization is applied:
    * Trial values of omegaI are sampled in the interval [omegaI_min, omegaI_max].
    * All MPI processes evaluate |W| for their assigned trial values.
    * The best candidate (minimum |W|) is selected at each iteration.
- The interval is repeatedly narrowed until delta_omegaI < delta_omegaI_tol.

## MPI Parallelization

The code is parallelized using MPI as follows:

- Each MPI process evaluates the Wronskian for a different trial value of omegaI.
- Values (omegaI, f(omegaI)) are gathered on rank 0.
- Rank 0 selects the best trial value:
    * If an improvement is found, the random-walk is restarted.
    * Otherwise, the search continues.
- Rank 0 broadcasts the updated omegaI to all processes.
- The loop continues until convergence.

This approach provides nearly linear speedup with the number of MPI processes.

## File Structure

bound_states_polar_fast.f90: Fast version of the code

bound_states_polar_exact.f90: Exact version of the code

README.txt: This file

## Requirements

- MPI Fortran compiler
- Quadruple precision support (REAL(16) and COMPLEX(16))

## Compilation

Compile using:

  mpiifort -O3 bound_states_polar.f90 -o bound_states_polar

or with optimizations:

  mpiifort -O3 -march=native -ffast-math bound_states_polar.f90 -o bound_states_polar

## Running the Program

Run with N MPI processes (Here, N is 8):

  mpirun -np 8 ./bound_states_polar

## Output Description

omegaI.txt:
    Contains the final converged value of omegaI and the absolute value of the Wronskian for it.

rZmZpW.txt:
    Contains r, Zm, Zp, W, and |W|.

## Notes

- Ensure that MPI is correctly installed.
- Quadruple precision is required for numerical stability.
- The search interval for omegaI may need adjustment for different l.

## REPOSITORY:

https://github.com/krezazadeh/bound_states_polar

## Reference

https://www.arxiv.org/abs/????.?????

## Contact

Kazem Rezazadeh

School of Astronomy,
Institute for Research in Fundamental Sciences (IPM)

Email: kazem.rezazadeh@ipm.ir

GitHub: https://github.com/krezazadeh

### Date

14 November 2025


