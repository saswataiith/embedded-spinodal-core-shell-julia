# Embedded spinodal decomposition: core-shell and Janus studies

These codes implement the embedded-particle model used in my research on phase separation and core-shell/Janus morphologies. The codes and the embedded-model idea are mine. I recovered the CUDA copies from my email correspondence with Pankaj; their location in that correspondence does not identify their author. I retain those copies in `original/`. The Julia implementation follows the 2022 equations. It uses the CPU on macOS, Linux or Windows and does not require a GPU. GPU acceleration is not yet implemented.

This readable CPU implementation follows the fixed-particle model in the recovered 13 July 2022 CUDA source. It is a new implementation of that source's equations, not a port of the different 2016 two-composition model or the evolving-particle 2023 implementation. Recovered source files in `original/` are unchanged. Historical executables and output logs are not included.

## Run

Julia 1.12.6 and FFTW were used for these checks. In this directory:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. core_shell.jl
julia --project=. checks.jl
```

The first command installs FFTW. The simulation writes `results-256/composition.csv`, `particle.csv` and `history.csv`. Rows of history contain time, mean composition, total free energy, minimum composition, maximum composition and instantaneous continuum-rate discrete dissipation. Fields are square arrays. Default run: 256² points, spacing 1, timestep 0.01, 1000 steps, final time 10. This is a short verification run; the recovered input requests 40 million steps. This simple CPU version allocates temporary FFT arrays and is intended for inspection before optimization.

## Equations and quantities

All parameters below are dimensionless; the recovered input does not supply a physical length, energy or time conversion. No SI units are inferred.

Composition c is the local composition variable. Particle field φ equals approximately 1 inside and 0 outside and is fixed throughout this version. Define h(φ)=φ³(10−15φ+6φ²) and g(φ)=φ²(1−φ)².

Bulk free-energy density:

f(c,φ)=A[1−h(φ)](c−cα)²+B h(φ)(c−cβ)²(c−cγ)²+P(1−χc)g(φ).

Total free energy F=∫[f+κ|∇c|²] dA. Therefore μ=∂f/∂c−2κ∇²c. The factor 2 follows from this gradient-energy convention. Mobility M=M₀h(φ). Composition evolves as ∂c/∂t=∇·(M∇μ), with periodic boundaries. There is no elasticity in this implementation.

A=2 sets the exterior well; B=4 sets the particle double well. Exterior preferred composition cα=0.5; particle wells cβ=0 and cγ=1. P=5 sets the surface barrier contribution, χ=0.5 sets its composition dependence, κ=2 sets gradient energy, M₀=1 sets mobility and Q=0.5 is the numerical stabilization parameter. Increasing c lowers the surface term for positive χ; this produces surface preference for the high-c phase. It does not by itself determine the final morphology.

Initial φ=[1−tanh((r−R)/w)]/2, R=70, w=4. Inside φ≥0.5, c=0.5 plus uniform noise of amplitude 0.01 with its interior mean removed; exterior c=0.5. Seed=5749. Julia's Mersenne Twister is reproducible but differs from cuRAND, so the original seed does not yield the same individual random field.

The Fourier update is ĉ(new)=ĉ+dt·FFT[∇·(M∇μ)]/[1+2Qκk⁴dt]. This matches the recovered stabilized update. The extra old-field term in its original numerator is intentional; it is not an erroneous cancellation. No unconditional energy-stability claim is made. First-derivative Nyquist modes are zero on even grids, as in the original; the energy and chemical potential retain the Laplacian Nyquist contribution.

## Findings and suggested corrections

1. The 2022 chemical-potential derivative is consistent with the free energy above. Its φ is fixed: relax_coeff and κφ do not evolve it. The attachment's phrase “for elasticity addition” should not be interpreted as an existing elastic solver.
2. Initial composition is written as c₀(cβ+cγ) in the original. This only implements the intended interpolation for the supplied wells 0 and 1. For a fractional position between arbitrary wells use (1−c₀)cβ+c₀cγ.
3. Original circular geometry uses grid indices without dx or dy, while another shape uses physical distances. Set a single length convention before any grid-refinement study. The Julia version explicitly uses physical nondimensional coordinates.
4. Two initialization/restart GPU kernels use tid≤nx·ny instead of tid<nx·ny. Change this to avoid an out-of-range element on some launch configurations. It is not necessarily triggered by the original exactly tiled launch.
5. The loop runs count=0 through total_steps inclusive and advances at the final count. Stop before the final update to make the requested final time and saved final field agree.
6. The active 2022 mobility ignores mob_tol; later noise injection is commented out. Label these settings inactive or implement them explicitly.
7. The original composition check only rejects c<−0.4 or c>1.4. Polynomial Cahn–Hilliard models do not guarantee c stays in [0,1]. Report bounds and investigate overshoots; do not clip c, because clipping changes conservation and the model.
8. A small change per timestep alone does not prove equilibrium. Include chemical-potential-gradient/dissipation diagnostics and repeat at a smaller dt. Degenerate mobility prevents transport in the far exterior; constant μ across the entire box is not the appropriate criterion there.
9. The 2023 3D source includes both conserved c and Allen–Cahn φ, with the extra exterior energy offset A₀. Its local derivatives are consistent with adding (1−h)A₀ to f. It mixes a spectral chemical-potential Laplacian with fourth-order finite-difference gradients/divergence. Such a combination is not automatically invalid: paired periodic first derivatives can retain the dissipation identity, but the identity and timestep behavior must be tested for that implementation.
10. Definite 2023 error: ComputeGradmu_z is passed dy in evolve.cu while the z divergence uses dz. Pass dz to both. Equal spacings conceal the error. The 2023 φ bounds check is disabled, and the coupled-field stopping condition requires further review. The original CUDA programs have not been compiled or run on this Mac; this report does not certify all their execution paths.

## Verification completed

`check-results.txt` records actual runs. Analytical μ agreed with a centered finite variation of F. The periodic discrete identity ∑μ·(∂c/∂t)=−∑M|∇μ|² passed. Mean composition remained constant to <1.3×10⁻¹⁵ in the tested runs, and energy decreased at every recorded step.

At 64², radius 16, final time 10, dt=0.01, 0.005 and 0.0025 gave successive field RMS differences 2.02746×10⁻⁴ and 1.02310×10⁻⁴, consistent with first-order time refinement for this case. This is not spatial convergence.

At 256², radius 70, final time 10, energy fell from 4406.510455 to 3592.020751, mean drift was 6.66×10⁻¹⁶, and final composition ranged from 0.00736 to 0.98997. The saved composition is a transient field, not an established equilibrium core-shell or Janus result.

Next: longer runs, equal physical particle and interface sizes across grid refinements, comparisons at the same physical time, and morphology measures at matched model parameters. A corrected coupled c–φ 3D port should be treated separately rather than silently substituted for this fixed-particle model.


## Repository contents

- `core_shell.jl`: readable Julia simulation.
- `checks.jl`: numerical checks.
- `Project.toml` and `Manifest.toml`: dependency files.
- `check-results.txt`: results from the completed checks.
- `original/2022-2d` and `original/2023-3d`: recovered CUDA sources and input files. These are references; the original 3D source still contains the errors described above.

The code and model attribution here is to my work. Published papers retain their original author lists. No new open-source licence is applied without a separate decision.


## Earlier related work

I used confined-system and embedded-domain approaches in my earlier research. Related publications include:

- *Phase separating bulk metallic glass: a hierarchical composite*, Physical Review Letters **96**, 245503 (2006), [doi:10.1103/PhysRevLett.96.245503](https://doi.org/10.1103/PhysRevLett.96.245503). This paper concerns core-shell and hierarchical phase-separating structures.
- *Anisotropic Li intercalation in a LiₓFePO₄ nano-particle: a spectral smoothed boundary phase-field model*, Physical Chemistry Chemical Physics **18**, 9537–9543 (2016), [doi:10.1039/C6CP00267F](https://doi.org/10.1039/C6CP00267F). This paper uses a spectral smoothed-boundary model for a particle in its surrounding medium.

These references describe earlier related research. The present repository is not claimed to reproduce either paper with its supplied short example.
