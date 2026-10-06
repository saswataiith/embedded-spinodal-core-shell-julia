# Embedded spinodal decomposition: core-shell and Janus studies

## Development of the model

In my [PRL work (2006)](https://doi.org/10.1103/PhysRevLett.96.245503), I used composition $c_A$ as an auxiliary variable in a ternary model to study hierarchical phase-separating structures. In my [PCCP work (2016)](https://doi.org/10.1039/C6CP00267F), I introduced a separate particle-domain field $\phi$ and used it to impose boundary conditions.

We studied core-shell precipitation in Sandeep Sugathan’s thesis, *A Phase-Field Study of Elastic Stress Effects on Phase Separation in Ternary Alloy Systems* (IIT Hyderabad, 2019), and in our [MRS Advances paper (2019)](https://doi.org/10.1557/adv.2019.104). That model couples conserved compositions to an evolving order parameter, includes elasticity, and uses $1-h(\phi)$ and $h(\phi)$ to interpolate between matrix and precipitate free energies.

The later core-shell and Janus studies with Pankaj and Subhradeep Chatterjee are described in [Acta Materialia (2022)](https://doi.org/10.1016/j.actamat.2022.117933), [MSMSE (2023)](https://doi.org/10.1088/1361-651X/aca420), and Pankaj’s thesis. I developed the embedded-particle model and the original CPU code. Pankaj developed the GPU implementation from my CPU code.

See [model history and references](MODEL-HISTORY.md) for the formulations and their differences.

## Present Julia implementation

The Julia code solves the two-dimensional Cahn–Hilliard equation in a fixed particle with composition-dependent surface energy and periodic boundaries. It runs on the CPU on macOS, Linux and Windows; Julia GPU acceleration is not implemented.

## Model

Composition $c$ describes the local alloy composition. The fixed particle field $\phi$ is approximately 1 inside and 0 outside. Define $h(\phi)=\phi^3(10-15\phi+6\phi^2)$ and $g(\phi)=\phi^2(1-\phi)^2$.

```math
f(c,\phi)=A[1-h(\phi)](c-c_\alpha)^2+B h(\phi)(c-c_\beta)^2(c-c_\gamma)^2+P(1-\chi c)g(\phi).
```

```math
F=\int_\Omega [f(c,\phi)+\kappa|\nabla c|^2]\,\mathrm{d}A.
```

```math
\mu=\frac{\partial f}{\partial c}-2\kappa\nabla^2 c,\qquad
\frac{\partial c}{\partial t}=\nabla\cdot[M_0 h(\phi)\nabla\mu].
```

Here $F$ is total free energy, $\mu$ is chemical potential, $\kappa$ is the gradient-energy coefficient and $M_0$ sets mobility. The coefficients $A$ and $B$ set the exterior and particle wells; $c_\alpha,c_\beta,c_\gamma$ are their preferred compositions. The parameters $P$ and $\chi$ control the surface-energy contribution. All quantities use nondimensional units. The solver uses a stabilized Fourier update; it does not include elasticity.

## Run

Install Julia, then run:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. core_shell.jl
julia --project=. checks.jl
```

The default example uses $256^2$ points, spacing 1, timestep 0.01 and 1000 steps (final time 10). Parameters and initialization are defined in `core_shell.jl`. Outputs in `results-256/` are `composition.csv`, `particle.csv` and `history.csv`; the history records time, mean composition, energy, composition bounds and dissipation.

## Checks and source versions

Julia 1.12.6 and FFTW checks covered the free-energy derivative, discrete dissipation identity, mass conservation and timestep refinement. In the $256^2$ run, mean-composition drift was $6.66\times10^{-16}$ and energy decreased from 4406.51 to 3592.02. This short run gives a transient field. Long-time morphology and spatial convergence remain to be tested. See [check results](check-results.txt) and [numerical details](docs/NUMERICAL-REVIEW.md).

Unchanged historical sources are in `original/`: `cpu-2019` uses constant mobility and evolves the particle field; `2022-2d` is the fixed-particle CUDA model implemented here; `2023-3d` evolves both fields. These historical programs have not been run in these checks. Checksums are in `source-manifest.json`.
