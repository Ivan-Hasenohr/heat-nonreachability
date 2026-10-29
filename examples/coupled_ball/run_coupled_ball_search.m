% Proposition 22: search for a witness separating an L2 ball.
% Initial state: (sin(pi*x),0).
% Target: the product-L2 ball centered at zero, with radius epsilon.
% target_radius: required scalar intval for a ball; ignored otherwise.
% Final time: T = 0.32.

project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root, 'src')));

%% Initialize the coarse discretization

Nx = 100;
Nt = 2000;

% Check the search-mesh limits before allocating the propagation table.
assert(Nx <= 100 && Nt <= 2000, ...
    'Use a coarse mesh for witness search.');

[args,args_float] = init_coupled_heat_fem( ...
    Nx,Nt,intval('0.32'),true);

expected_T = intval('0.32');
assert(inf(args{1}) == inf(expected_T) ...
    && sup(args{1}) == sup(expected_T), ...
    'The final time must be exactly intval(''0.32'').');


%% Define the exact continuous states and target set

% Sine coefficients, concatenated component by component.
y_0_base = intval([1; 0]);
y_f_base = intval([0; 0]);

% Both vectors contain concatenated sine-series coefficients.
n_coeffs = [length(y_0_base), length(y_f_base)];
assert(all(n_coeffs ~= 2*(Nx-1)) ...
    && all(n_coeffs ~= 4*Nx), ...
    'The sine data would be interpreted as P1 or Hermite data on this mesh.');

epsilon = intval(1)/100;


target_kind = 2;

% Floating samples used only for witness search.
x = (1:Nx-1)'/Nx;
y_0 = [sin(pi*x); zeros(Nx-1,1)];
y_f = zeros(2*Nx-2,1);

%% Initialize the search

% Start from the target minus the uncontrolled final state.
y_tilde_T = y_f - args_float{8}(:,:,Nt+1)*y_0;

p_0 = y_tilde_T;
assert(norm(p_0) > 0, ...
    'The initial witness must be nonzero.');

p_0 = p_0/(norm(p_0)*sqrt(1/Nx));
u_0 = zeros(Nx-1,Nt);

%% Search and evaluate the resulting witness

started = tic;

[J, p, u] = find_and_certify_coupled_heat( ...
    y_0, y_0_base, y_f, y_f_base, ...
    p_0, u_0, args, args_float, target_kind, epsilon);

search_and_certification_seconds = toc(started);

%% Preserve the witness independently of certification success

certified = isfinite(sup(J)) && sup(J) < 0;

save(fullfile(project_root, 'data', 'witnesses', ...
    'proposition22_coarse.mat'), ...
    'p', 'u', 'Nx', 'Nt', 'args', ...
    'y_0_base', 'y_f_base', 'epsilon', 'target_kind', ...
    'J', 'certified', 'search_and_certification_seconds');

%% Display floating search diagnostics

control_field = fl_apply_control_adjoint(p, args_float);


% Floating diagnostic: integrated P1 control support with nodal target terms.
% Certification is reported separately below.
J_search = fl_p1_box_field_support(control_field, args_float) ...
    + fl_ball_target_term(-p, y_tilde_T, mid(epsilon));

disp('Final floating diagnostic (P1 control support, nodal target terms):');
disp(J_search);

disp('Coarse certification enclosure:');
disp(J);

disp('Certified non-reachability on this mesh:');
disp(certified);