# Projected scalar heat target

`run_heat_projected_target_certificate.m` certifies Proposition 22 with
T=13/100, zero initial state and target P6 L_T(3/4). It encloses the exact
sine coefficients by INTLAB and uses the fixed 24-mode input
`data/witnesses/heat_projected_sine.mat`.

The mesh is Nx=2000, Nt=1200000. The verified witness is the natural spline
of the fixed normalized nodal samples. The result is saved to
`results/heat_projected.mat`. Run the complete campaign to obtain this
result alongside all other examples from one frozen source tree.
