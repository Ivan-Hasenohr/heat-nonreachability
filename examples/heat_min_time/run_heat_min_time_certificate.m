% Corollary 19: certify the witness optimized specifically for T = 0.13.

project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root, 'src')));

% Create the output directory on a fresh checkout.
if ~isfolder(fullfile(project_root,'results'))
    mkdir(fullfile(project_root,'results'));
end

%% Initialize the fine discretization
args = init_heat_fem(1200,400000,intval('0.13'),false);
Nx = args{11};

%% Transfer the newly optimized coarse witness

coarse = load(fullfile(project_root, 'data', 'witnesses', ...
    'corollary19_coarse.mat'));

p = resample_nodal_witness(coarse.p, Nx);

% Normalize using the nodal norm.
% This does not impose an exact unit L2 norm on the continuous spline.
p = p/(norm(p)*sqrt(1/Nx));

%% Define the exact continuous target

y_f_base = intval(1)/50;

% This target is a sine series, not P1 or cubic Hermite data.
assert(length(y_f_base) ~= Nx-1 ...
    && length(y_f_base) ~= 2*Nx, ...
    'The sine target would be interpreted as P1 or Hermite data on this mesh.');

y_0 = zeros(Nx-1,1);

% Preserve the fine witness before certification.
save(fullfile(project_root, 'data', 'witnesses', ...
    'corollary19_fine.mat'), 'p', 'Nx', 'y_f_base');

%% Certify and save the result

started = tic;

[certified, J, support, pairing, err, bounds] = ...
    certify_heat_nonreachability(p, y_0, y_f_base, args);

certification_seconds = toc(started);

save(fullfile(project_root, 'results', 'corollary19.mat'), ...
    'p', 'y_0', 'y_f_base', 'args', ...
    'certified', 'J', 'support', 'pairing', 'err', 'bounds', 'certification_seconds');

disp('Certified non-reachability:');
disp(certified);

disp('Final certification enclosure:');
disp(J);