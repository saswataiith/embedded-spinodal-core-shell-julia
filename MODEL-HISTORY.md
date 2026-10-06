# How I developed the auxiliary-variable models

In my PRL work, I used the composition cA as the auxiliary variable in my ternary model. The three composition variables obey cA+cB+cC=1. See *Phase separating bulk metallic glass: a hierarchical composite*, Physical Review Letters 96, 245503 (2006), https://doi.org/10.1103/PhysRevLett.96.245503.

In my PCCP work, I introduced a separate auxiliary field phi and used it to impose the boundary conditions. See *Anisotropic Li intercalation in a LiₓFePO₄ nano-particle: a spectral smoothed boundary phase-field model*, Physical Chemistry Chemical Physics 18, 9537–9543 (2016), https://doi.org/10.1039/C6CP00267F.

In the present embedded-particle model, I use 1−h(phi) and h(phi) to weight the exterior and particle free-energy contributions. Here h(phi)=phi³(10−15phi+6phi²). Composition and phi are dimensionless. The field phi distinguishes the two regions; h smoothly interpolates between their free energies.

We explored this interpolation in Sandeep Sugathan’s thesis, *A Phase-Field Study of Elastic Stress Effects on Phase Separation in Ternary Alloy Systems* (IIT Hyderabad, June 2019). Chapter 3, printed page 45, Eq. (3.77), gives:

f0=(1−h(phi))f_alpha+h(phi)f_beta_gamma+Omega phi²(1−phi)².

Here f_alpha is the matrix free energy, f_beta_gamma describes the two ordered precipitate phases, and Omega sets the barrier between matrix and precipitate. Omega has the same free-energy scale as f0. Eq. (3.76) includes composition and phi gradient energies.

Chapter 6, Section 6.1, printed page 155, applies this formulation to core-shell precipitation. Phi=0 describes the disordered alpha matrix; phi=1 describes either ordered precipitate, beta or gamma. Their compositions distinguish the two precipitate phases. In that thesis model, phi evolves by an Allen–Cahn equation (Chapter 3, Eq. (3.80)). In the 2022 embedded-particle code, phi remains fixed.

These models use different free energies and evolution equations. I developed the embedded-model idea and the original CPU code. Pankaj developed the GPU implementation from my CPU code. Our core-shell and Janus studies with Subhradeep Chatterjee are described in [Acta Materialia (2022)](https://doi.org/10.1016/j.actamat.2022.117933), [Modelling and Simulation in Materials Science and Engineering (2023)](https://doi.org/10.1088/1361-651X/aca420), and Pankaj’s doctoral thesis.
