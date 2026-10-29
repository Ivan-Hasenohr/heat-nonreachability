% Proposition 22: certify a previously computed two-component witness.
% Target: the product-L2 ball centered at zero, with radius epsilon.
% target_radius: required scalar intval for a ball; ignored otherwise.

project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root, 'src')));

% Create the output directory on a fresh checkout.
if ~isfolder(fullfile(project_root,'results'))
    mkdir(fullfile(project_root,'results'));
end

%% Initialize the fine discretization
Nx = 500;
Nt = 25000;
args = init_coupled_heat_fem(Nx,Nt,intval('0.32'),false);

%% Load the coarse witness

coarse = load(fullfile(project_root, 'data', 'witnesses', ...
    'proposition22_coarse.mat'));

coarse_Nx = coarse.Nx;

assert(length(coarse.p) == 2*(coarse_Nx-1), ...
    'The coarse witness must contain two concatenated components.');

%% Transfer each component to the fine mesh

p_1 = resample_nodal_witness( ...
    coarse.p(1:coarse_Nx-1), Nx);

p_2 = resample_nodal_witness( ...
    coarse.p(coarse_Nx:end), Nx);

% Normalize jointly using the nodal norm, preserving relative component amplitudes.
% The continuous spline L2 norm is evaluated separately during certification.
p = [p_1; p_2];

assert(norm(p) > 0, ...
    'The transferred witness must be nonzero.');

p = p/(norm(p)*sqrt(1/Nx));

%% Define the exact continuous states and target set

y_0_base = intval([1; 0]);
y_f_base = intval([0; 0]);

% Both vectors contain concatenated sine-series coefficients.
n_coeffs = [length(y_0_base), length(y_f_base)];
assert(all(n_coeffs ~= 2*(Nx-1)) ...
    && all(n_coeffs ~= 4*Nx), ...
    'The sine data would be interpreted as P1 or Hermite data on this mesh.');

epsilon = intval(1)/100;

target_kind = 2;

%% Preserve the fine witness before certification

save(fullfile(project_root, 'data', 'witnesses', ...
    'proposition22_fine.mat'), ...
    'p', 'Nx', 'Nt', 'y_0_base', 'y_f_base', ...
    'epsilon', 'target_kind');

%% Certify non-reachability of the entire target ball

started = tic;

[certified, J, bounds] = certify_coupled_heat_nonreachability( ...
    p, y_0_base, y_f_base, args, target_kind, epsilon);

certification_seconds = toc(started);

certified = logical(certified) ...
    && isfinite(sup(J)) && sup(J) < 0;

%% Save the result independently of certification success

save(fullfile(project_root, 'results', 'proposition22.mat'), ...
    'p', 'args', 'Nx', 'Nt', 'y_0_base', 'y_f_base', ...
    'epsilon', 'target_kind', 'certified', 'J', ...
    'certification_seconds', 'bounds');

disp('Certified non-reachability:');
disp(certified);

disp('Final certification enclosure:');
disp(J);

disp('Certification time in seconds:');
disp(certification_seconds);