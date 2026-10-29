% Proposition 20: regularized witness search for the tent target.

project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root, 'src')));

%% Initialize the coarse discretization


[args,args_float] = init_heat_fem(100,2000,intval('1'),true);

Nx = args{11};
Nt = args{2};

%% Define the exact continuous target

% The tent has its vertex at (4/5, 1/40).
% Its P1 representation is exact when Nx is divisible by 5.
assert(mod(Nx,5) == 0, 'Nx must be divisible by 5.');

j = (1:Nx-1)';
left = 5*j <= 4*Nx;

% y_f(x) = x/32 for x <= 4/5.
% y_f(x) = (1-x)/8 for x >= 4/5.
y_f_base = intval(zeros(Nx-1,1));

y_f_base(left) = intval(j(left))/(32*intval(Nx));

y_f_base(~left) = (intval(Nx)-intval(j(~left)))/(8*intval(Nx));


% This target is represented by P1 interior nodal values.
assert(length(y_f_base) == Nx-1, ...
    'The tent target must contain Nx-1 P1 nodal values.');


%% Initialize the floating search

y_f = mid(y_f_base);
p_0 = y_f/(norm(y_f)*sqrt(1/Nx));
u_0 = zeros(Nx-1,Nt);

lambda_reg = 1e-12;

%% Search and certify the resulting witness

% Regularization affects witness search only.
% The verifier evaluates the original, unregularized functional.
[J, p, u] = find_and_certify_heat( ...
    y_f, y_f_base, p_0, u_0, args, args_float, lambda_reg);

%% Preserve the witness even if coarse certification is inconclusive

save(fullfile(project_root, 'data', 'witnesses', 'proposition20_coarse.mat'), 'p', 'Nx', 'Nt', 'lambda_reg', 'y_f_base', 'J', 'args');

%% Display diagnostics

J_search = fl_heat_positive_control_support(p, args_float) - p'*y_f/Nx;

disp('Final floating unregularized functional:');
disp(J_search);

disp('Coarse certification enclosure:');
disp(J);