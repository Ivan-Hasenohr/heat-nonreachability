% Corollary 18: certify a previously saved witness on a fine mesh.

project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root, 'src')));

% Create the output directory on a fresh checkout.
if ~isfolder(fullfile(project_root,'results'))
    mkdir(fullfile(project_root,'results'));
end

% Load the witness obtained with both control intervals.
coarse = load(fullfile(project_root, 'data', 'witnesses', ...
    'corollary18_coarse.mat'));

args = init_heat_fem(1200,3200000,intval('1'),false);

Nx = args{11};

% Transfer the coarse natural cubic spline to the fine mesh.
p = resample_nodal_witness(coarse.p, Nx);

% Normalize using the nodal norm.
% This does not impose an exact unit L2 norm on the continuous spline.
p = p/(norm(p)*sqrt(1/Nx));

% Preserve the exact continuous target.
y_f_base = intval('0.035');

% This target is a sine series, not P1 or cubic Hermite data.
assert(length(y_f_base) ~= Nx-1 ...
    && length(y_f_base) ~= 2*Nx, ...
    'The sine target would be interpreted as P1 or Hermite data on this mesh.');


y_0 = zeros(Nx-1,1);

% Save the fine witness before the potentially long certification.
save(fullfile(project_root, 'data', 'witnesses', ...
    'corollary18_fine.mat'), 'p', 'Nx', 'y_f_base');

started = tic;


[certified,J,support,pairing,err,bounds] = ...
    certify_heat_nonreachability(p,y_0,y_f_base,args);

certification_seconds = toc(started);

save(fullfile(project_root, 'results', 'corollary18.mat'), ...
    'p', 'y_0', 'y_f_base', 'args', ...
    'certified', 'J', 'support', 'pairing', 'err', 'bounds', "certification_seconds");

disp('Final certification enclosure:');
disp(J);