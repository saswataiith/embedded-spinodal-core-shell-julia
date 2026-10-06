# Development of my auxiliary-variable models

## PRL: an auxiliary composition in the ternary model

In my PRL work, I used composition $c_A$ as the auxiliary variable in a ternary model. The dimensionless local compositions satisfy $c_A+c_B+c_C=1$. This work concerns hierarchical phase-separating structures.

**Reference:** B. J. Park et al., *Phase Separating Bulk Metallic Glass: A Hierarchical Composite*, Physical Review Letters **96**, 245503 (2006). [Paper](https://doi.org/10.1103/PhysRevLett.96.245503).

## PCCP: a separate particle-domain field

In my PCCP work, I introduced a separate auxiliary field $\phi$ and used it to impose the particle boundary conditions in a spectral smoothed-boundary model. This is distinct from using a ternary composition as the auxiliary variable.

**Reference:** L. Hong, L. Liang, S. Bhattacharyya, W. Xing and L. Q. Chen, *Anisotropic Li intercalation in a LiₓFePO₄ nano-particle: a spectral smoothed boundary phase-field model*, Physical Chemistry Chemical Physics **18**, 9537–9543 (2016). [Paper](https://doi.org/10.1039/C6CP00267F).

## Sandeep’s thesis and the MRS Advances core-shell study

We explored matrix–precipitate free-energy interpolation in Sandeep Sugathan’s thesis, *A Phase-Field Study of Elastic Stress Effects on Phase Separation in Ternary Alloy Systems* (IIT Hyderabad, June 2019). Chapter 3, Eq. (3.77), printed page 45, gives

```math
f_0=(1-h(\phi))f_\alpha+h(\phi)f_{\beta\gamma}+\Omega\phi^2(1-\phi)^2,
\qquad h(\phi)=\phi^3(10-15\phi+6\phi^2).
```

Here $f_\alpha$ is the matrix free-energy density; $f_{\beta\gamma}$ describes the two ordered precipitate phases; and $\Omega$ sets the barrier between matrix and precipitate. These three quantities have the same energy-density units before nondimensionalization. The order parameter $\phi$ and interpolation function $h$ are dimensionless. Equation (3.76) includes composition and order-parameter gradient energies.

Chapter 6, Section 6.1, printed page 155, applies this formulation to core-shell precipitation. The value $\phi=0$ identifies the disordered matrix $\alpha$; $\phi=1$ identifies either ordered precipitate phase, $\beta$ or $\gamma$. Their compositions distinguish the precipitate phases. Conserved compositions evolve through Cahn–Hilliard equations, while $\phi$ evolves through an Allen–Cahn equation. Elastic energy is included.

**Related paper:** S. Sugathan and S. Bhattacharya, *Phase-Field Modelling of Evolution of Compact Ordered Precipitates in Ternary Alloy Systems*, MRS Advances **4**(25–26), 1457–1463 (2019). [Paper](https://doi.org/10.1557/adv.2019.104).

## Later core-shell and Janus studies

I developed the embedded-particle model and the original CPU code. Pankaj developed the GPU implementation from my CPU code. Our studies with Subhradeep Chatterjee examine the competition between phase separation and surface effects in bimetallic particles.

- P. Pankaj, S. Bhattacharyya and S. Chatterjee, *Competition of core-shell and Janus morphology in bimetallic nanoparticles: Insights from a phase-field model*, Acta Materialia **233**, 117933 (2022). [Paper](https://doi.org/10.1016/j.actamat.2022.117933) · [Preprint](https://arxiv.org/abs/2105.09531).
- P. Pankaj, S. Bhattacharyya and S. Chatterjee, *Surface-directed and bulk spinodal decomposition compete to decide the morphology of bimetallic nanoparticles*, Modelling and Simulation in Materials Science and Engineering **31**, 015003 (2023). [Paper](https://doi.org/10.1088/1361-651X/aca420) · [Preprint](https://arxiv.org/abs/2206.08274).
- Pankaj’s doctoral thesis, IIT Hyderabad: embedded-domain formulation and nanoparticle morphology studies.

## Scope of this Julia code

The current Julia implementation follows the fixed-particle 2022 model. It uses $1-h(\phi)$ and $h(\phi)$ to weight the exterior and particle free energies, but keeps $\phi$ fixed and evolves only composition. It does not include elasticity. The historical 2019 CPU and 2023 CUDA implementations evolve the particle field and are retained separately. These models share the auxiliary-variable development described above but have different free energies, kinetic equations and boundary treatments. The short Julia example is not a reproduction of the earlier publications.
