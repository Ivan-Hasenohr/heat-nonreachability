% Proposition 20: certify a previously computed regularized witness.
% Certification evaluates the original, unregularized functional.

project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root, 'src')));

% Create the output directory on a fresh checkout.
if ~isfolder(fullfile(project_root,'results'))
    mkdir(fullfile(project_root,'results'));
end

%% Initialize the fine discretization
Nx = 1200;
Nt = 3200000;
args = init_heat_fem(Nx,Nt,intval('1'),false);

%% Load and transfer the regularized witness

coarse = load(fullfile(project_root, 'data', 'witnesses', ...
    'proposition20_coarse.mat'));

lambda_reg = coarse.lambda_reg;

p = resample_nodal_witness(coarse.p, Nx);

% Normalize using the nodal norm.
% This does not impose an exact unit L2 norm on the continuous spline.
p = p/(norm(p)*sqrt(1/Nx));

%% Reconstruct the exact continuous tent target

% The vertex (4/5,1/40) must be a mesh node.
assert(mod(Nx,5) == 0, ...
    'Nx must be divisible by 5.');

j = (1:Nx-1)';
left = 5*j <= 4*Nx;

% y_f(x) = x/32 for x <= 4/5.
% y_f(x) = (1-x)/8 for x >= 4/5.
y_f_base = intval(zeros(Nx-1,1));

y_f_base(left) = ...
    intval(j(left))/(32*intval(Nx));

y_f_base(~left) = ...
    (intval(Nx)-intval(j(~left)))/(8*intval(Nx));

% This target is represented by P1 interior nodal values.
assert(length(y_f_base) == Nx-1, ...
    'The tent target must contain Nx-1 P1 nodal values.');


% Floating samples are not the rigorous definition of the target.
y_0 = zeros(Nx-1,1);

%% Save the fine witness before certification

save(fullfile(project_root, 'data', 'witnesses', 'proposition20_fine.mat'), 'p', 'Nx', 'Nt', 'lambda_reg', 'y_f_base');

%% Certify the original functional

started = tic;

[certified,J,support,pairing,err,bounds] = ...
    certify_heat_nonreachability(p,y_0,y_f_base,args);

certification_seconds = toc(started);

% Require a finite, strictly negative upper bound.
certified = logical(certified) ...
    && isfinite(sup(J)) && sup(J) < 0;

%% Save the full result

save(fullfile(project_root, 'results', 'proposition20.mat'), 'p', 'y_0', 'y_f_base', 'args', 'lambda_reg', 'certified', 'J', 'support', 'pairing', 'err', 'certification_seconds', "bounds");

disp('Certified non-reachability:');
disp(certified);

disp('Final certification enclosure:');
disp(J);

disp('Certification time in seconds:');
disp(certification_seconds);