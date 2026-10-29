% Corollary 19: search for a witness at T = 0.13.
% Start from the coarse witness obtained for Corollary 18.

project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root, 'src')));

%% Initialize the coarse discretization

[args,args_float] = init_heat_fem(100,2000,intval('0.13'),true);
Nx = args{11};
Nt = args{2};

%% Define the target

% Exact continuous target: sin(pi*x)/50.
y_f_base = intval(1)/50;

% This target is a sine series, not P1 or cubic Hermite data.
assert(length(y_f_base) ~= Nx-1 ...
    && length(y_f_base) ~= 2*Nx, ...
    'The sine target would be interpreted as P1 or Hermite data on this mesh.');

% Floating samples used only for witness search.
x = (1:Nx-1)'/Nx;
y_f = mid(y_f_base)*sin(pi*x);

%% Initialize the search

% Reuse the coarse witness from Corollary 18.
previous = load(fullfile(project_root, 'data', 'witnesses', ...
    'corollary18_coarse.mat'));

assert(length(previous.p) == Nx-1, ...
    'The initial witness must match the search mesh.');

p_0 = previous.p(:);
p_0 = p_0/(norm(p_0)*sqrt(1/Nx));

u_0 = zeros(Nx-1,Nt);

%% Search and evaluate the resulting witness

[J, p, u] = find_and_certify_heat( ...
    y_f, y_f_base, p_0, u_0, args, args_float);

%% Save the witness independently of certification success

save(fullfile(project_root, 'data', 'witnesses', ...
    'corollary19_coarse.mat'), ...
    'p', 'Nx', 'Nt', 'y_f_base', 'J', 'args');

%% Display diagnostics

J_search = fl_heat_positive_control_support(p, args_float) ...
    - p'*y_f/Nx;

disp('Final floating search functional:');
disp(J_search);

disp('Coarse certification enclosure:');
disp(J);