% Corollary 18: search for a dual witness on a coarse mesh.

project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root, 'src')));

[args,args_float] = init_heat_fem(100,2000,intval('1'),true);

Nx = args{11};
Nt = args{2};

% Exact continuous target: 0.035*sin(pi*x).
y_f_base = intval('0.035');

% This target is a sine series, not P1 or cubic Hermite data.
assert(length(y_f_base) ~= Nx-1 ...
    && length(y_f_base) ~= 2*Nx, ...
    'The sine target would be interpreted as P1 or Hermite data on this mesh.');

% Floating target samples and initial guesses.
x = (1:Nx-1)'/Nx;
y_f = mid(y_f_base)*sin(pi*x);
p_0 = y_f/(norm(y_f)*sqrt(1/Nx));
u_0 = zeros(Nx-1,Nt);

[J, p, u] = find_and_certify_heat(y_f, y_f_base, p_0, u_0, args, args_float);

% Preserve the witness even if coarse certification is inconclusive.
save(fullfile(project_root, 'data', 'witnesses','corollary18_coarse.mat'), 'p', 'Nx', 'Nt', 'y_f_base', 'J');

disp('Coarse certification enclosure:');
disp(J);