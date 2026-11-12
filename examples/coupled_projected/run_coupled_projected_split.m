% Proposition 27: P4 L_T u_ref, where u_ref is 1 and -1/2 on the two subintervals.
% INTLAB must already be initialized.
example_dir = fileparts(mfilename('fullpath'));
addpath(example_dir);

Nx = 3090;
Nt = 47000;

result = run_coupled_projected_certificate('split',Nx,Nt,'fine');
disp(result);
