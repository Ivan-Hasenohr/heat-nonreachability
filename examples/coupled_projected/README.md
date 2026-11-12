# Coupled heat examples

After INTLAB initialization, add this directory to the MATLAB path.

| Script | Paper result | Target |
|---|---|---|
| run_coupled_projected_ball.m | Proposition 24 | Ball centred at P2 L_T(9/5), radius 1/5000, T=8/25 |
| run_coupled_projected_partial.m | Proposition 26 | {-23 sin(pi x)/10000 + sin(2 pi x)/1000} x L2, T=1/5 |
| run_coupled_projected_split.m | Proposition 27 | P4 L_T(u_ref), u_ref=1 on (0,1/4) and -1/2 on (1/4,1/2), T=3/50 |

All cases start at zero and use the same coupled diffusion/reaction operator,
control region (0,1/2) and control bounds [-1,2]. The shared function
`run_coupled_projected_certificate(name,Nx,Nt,run_name)` defines the exact
data, samples the fixed seed, constructs verified FEM arguments and calls
the existing interval certifier. The split target is evaluated directly at
its two control amplitudes. Its fixed seed is `data/witnesses/coupled_split_sine.csv`.

Outputs are saved under `results/coupled_projected/<name>/<run_name>.mat`
with a JSON summary. The same scripts run inside the frozen campaign tree.
