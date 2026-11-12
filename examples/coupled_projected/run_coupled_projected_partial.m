% First target component fixed; the second component is free, T=1/5.
% INTLAB must already be initialized. Edit Nx and Nt for a different mesh.
example_dir = fileparts(mfilename('fullpath'));
addpath(example_dir);

Nx = 1490;
Nt = 37000;

result = run_coupled_projected_certificate('partial',Nx,Nt,'fine');
disp(result);
