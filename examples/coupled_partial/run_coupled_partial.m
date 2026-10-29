% Proposition 24: certify non-reachability of the partial target.

project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root, 'src')));

% Create the output directory on a fresh checkout.
if ~isfolder(fullfile(project_root,'results'))
    mkdir(fullfile(project_root,'results'));
end

args = init_coupled_heat_fem(1000,80000,intval('1'),false);
Nx = args{11};

% Exact sine coefficients, concatenated component by component.
% Initial state: (0,0).
% Target: (-sin(2*pi*x)/50,0), with the second component free.
y_0_base = intval([0; 0]);
y_f_base = [intval(0); -intval(1)/50; intval(0); intval(0)];

% Both vectors contain concatenated sine-series coefficients.
n_coeffs = [length(y_0_base), length(y_f_base)];
assert(all(n_coeffs ~= 2*(Nx-1)) ...
    && all(n_coeffs ~= 4*Nx), ...
    'The sine data would be interpreted as P1 or Hermite data on this mesh.');

% Floating coefficients define the fixed nodal witness.
% Its second component must remain exactly zero.
x = (1:Nx-1)'/Nx;
p = [-sin(2*pi*x); zeros(Nx-1,1)];

% Normalize using the nodal norm.
% This does not impose an exact unit L2 norm on the continuous spline.
p = p/(norm(p)*sqrt(1/Nx));

started = tic;

[certified, J, bounds] = certify_coupled_heat_nonreachability(p, y_0_base, y_f_base, args);

certification_seconds = toc(started);

save(fullfile(project_root,'results','proposition24.mat'), ...
    'p','args','y_0_base','y_f_base','certified','J','bounds', 'certification_seconds');

disp('Certified non-reachability:');
disp(certified);
disp('Final enclosure:');
disp(J);