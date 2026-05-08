The Stability of Prendergast Magnetic Fields - Outline
<br> Andrew Koren

1. Introduction
   - Magnetic fields in stars: two origin hypotheses (dynamo vs. fossil), importance for massive stars and red giant cores.
   - The Prendergast field: an axisymmetric mix of poloidal and toroidal components proposed as a stable fossil field configuration.
   - Previous work (including Prendergast 1958) suggested short-term stability
   - We show Prendergast is unstable regardless of conditions / parameters

2. Numerical Details
   - Linearized MHD equations in Boussinesq approximation, spherical geometry, using Dedalus pseudospectral code.
   - Length by stellar core radius $$R$$, time by Alfvén time, magnetic field by max $$|B_0|$$.
   - Background field: Prendergast magnetic field. Using $$\lambda \approx 5.76346$$ as the lowest energy solution that is non-singular, vanishes at outer boundary.
   - Initial perturbation: azimuthal wavenumber $$m$$ selected via prescribed flow followed by divergence cleaning.
   - Boundary conditions: potential (POT) and perfectly conducting (PC); stress-free, impenetrable, no density perturbation at outer boundary.
   - Visualize magnetic field and velocity flow using Visualization and Analysis Platform for Ocean, Atmosphere, and Solar REsearchers (VAPOR)s

3. Results
    - The Prendergast magnetic field exhibits a robust linear instability with exponentially growing kinetic energy, independent of resolution and timestep, and its unstable eigenmode has an $$m=1$$ symmetry distinct from the background field.
    - The instability is resistive, with a growth rate that follows a power-law decrease as the magnetic resistivity is reduced.
    - Adding stable stratification reduces radial motion but does not stabilize the mode; the growth rate becomes constant at stratification strengths far below stellar values and the magnetic field perturbation remains unchanged.
    - Only azimuthal wavenumbers $$m=0$$ and $$m=1$$ are unstable, and the instability persists under both potential and perfectly conducting boundary conditions, though the growth rates and their power-law scalings differ.
    - The $$m=1$$ mode shows tearing-mode-like behavior at the origin, but tearing is not the sole driver of instability because the $$m=0$$ mode is also unstable without such a feature.
    - Extrapolating the resistive power-law scaling to stellar conditions yields an e-folding time of roughly $$10^6$$ years, shorter than the red giant branch lifetime, making the Prendergast field unsuitable as a long-lived magnetic equilibrium.

4. Conclusion
    - The Prendergast field is physically and robustly unstable due to a resistive mechanism that operates even at stellar stratifications and is fastest for the $$m=1$$ mode.
    - Because the instability grows on a timescale much shorter than stellar evolution, the Prendergast configuration cannot serve as a realistic model for fossil magnetic fields in stars.
    - Future work must investigate the nonlinear saturation of this instability and search for alternative truly stable magnetic equilibria.