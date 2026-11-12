% Ball centered at P2 L_T(9/5), radius 1/5000, T=8/25.
% INTLAB must already be initialized. Edit Nx and Nt for a different mesh.
example_dir = fileparts(mfilename('fullpath'));
addpath(example_dir);

Nx = 2220;
Nt = 129000;

result = run_coupled_projected_certificate('ball',Nx,Nt,'fine');
disp(result);
