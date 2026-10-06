# Embedded spinodal decomposition: core-shell and Janus studies

I developed the embedded-particle model and the original CPU code for studying phase separation in alloy particles. Pankaj developed the GPU implementation from my CPU code. Our related papers with Subhradeep Chatterjee are listed below.

The Julia code solves the two-dimensional Cahn–Hilliard equation in a fixed particle with composition-dependent surface energy and periodic boundaries. It runs on the CPU on macOS, Linux and Windows; Julia GPU acceleration is not implemented.

## Model

Composition $c$ describes the local alloy composition. The fixed particle field $\phi$ is approximately 1 inside and 0 outside. Define $h(\phi)=\phi^3(10-15\phi+6\phi^2)$ and $g(\phi)=\phi^2(1-\phi)^2$.

$$
f=A[1-h(\phi)](c-c_\alpha)^2+B h(\phi)(c-c_\beta)^2(c-c_\gamma)^2+P(1-\chi c)g(\phi).
$$

$$
F=\int_\Omega [f+\kappa|\nabla c|^2],dA,\qquad
\mu=\frac{\partial f}{\partial c}-2\kappa\nabla^2c,\qquad
\frac{\partial c}{\partial t}=\nabla\cdot[M_0h(\phi)\nabla\mu].
$$

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

## References

- P. Pankaj, S. Bhattacharyya and S. Chatterjee, **Competition of core-shell and Janus morphology in bimetallic nanoparticles: Insights from a phase-field model**, *Acta Materialia* **233**, 117933 (2022). [Paper](https://doi.org/10.1016/j.actamat.2022.117933) · [Preprint](https://arxiv.org/abs/2105.09531).
- P. Pankaj, S. Bhattacharyya and S. Chatterjee, **Surface-directed and bulk spinodal decomposition compete to decide the morphology of bimetallic nanoparticles**, *Modelling and Simulation in Materials Science and Engineering* **31**, 015003 (2023). [Paper](https://doi.org/10.1088/1361-651X/aca420) · [Preprint](https://arxiv.org/abs/2206.08274).
- Pankaj’s doctoral thesis, IIT Hyderabad: embedded-domain formulation and nanoparticle morphology studies.

For the earlier PRL and PCCP models and the formulation in Sandeep’s thesis, see [auxiliary-variable models](MODEL-HISTORY.md).
