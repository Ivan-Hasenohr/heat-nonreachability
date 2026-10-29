% Theorem 6: certify non-reachability at T = 1.15.

project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root,'src')));

%% Choose the certification mesh

Nx = 800;
Nt = 500000;
intval_T = intval('1.15');

validateattributes(Nx,{'double'},{'scalar','real','finite','integer','>=',3});
validateattributes(Nt,{'double'},{'scalar','real','finite','integer','>=',2});

fprintf('Certification mesh: Nx = %d, Nt = %d\n',Nx,Nt);
preparation_started = tic;

%% Transfer the saved witness to the chosen spatial mesh

source_witness = fullfile('data','witnesses','corollary18_fine.mat');
witness = load(fullfile(project_root,source_witness),'p','Nx');

assert(iscolumn(witness.p) && length(witness.p) == witness.Nx-1, ...
    'The saved witness dimensions do not match its spatial mesh.');

source_Nx = witness.Nx;
fprintf('Transferring the witness: Nx = %d -> %d...\n',source_Nx,Nx);

if Nx == source_Nx
    p = witness.p;
else
    p = resample_nodal_witness(witness.p,Nx);

    % Normalize using the nodal norm.
    % This does not impose an exact unit L2 norm on the continuous spline.
    p = p/(norm(p)*sqrt(1/Nx));
end

assert(iscolumn(p) && length(p) == Nx-1 && all(isfinite(p)) ...
    && norm(p) > 0,'The transferred witness is invalid.');

%% Define the verified problem

fprintf('Initializing the verified problem...\n');
args = init_heat_fem(Nx,Nt,intval_T,false);

y_f_base = intval('0.037');

% This target is a sine series, not P1 or cubic Hermite data.
assert(length(y_f_base) ~= Nx-1 ...
    && length(y_f_base) ~= 2*Nx, ...
    'The sine target would be interpreted as P1 or Hermite data on this mesh.');

y_0 = zeros(Nx-1,1);

%% Display the discretization budget before the long computation

p_base = iv_natural_cubic_from_nodes(p);
norm_y_f = norm(y_f_base)/intval('sqrt2');
discretization_bound = iv_heat_discretization_bound( ...
    p_base,y_0,norm_y_f,args,1,true);

preparation_seconds = toc(preparation_started);
fprintf('Total discretization bound: %.8g\n',discretization_bound);
fprintf('Witness preparation and initialization: %.1f s\n',preparation_seconds);

%% Certify the fixed witness

fprintf('Starting certification on the chosen mesh...\n');
started = tic;

[certified,J,support,pairing,err,bounds] = ...
    certify_heat_nonreachability(p,y_0,y_f_base,args);

certification_seconds = toc(started);
certified = isfinite(inf(J)) && isfinite(sup(J)) && certified;

%% Save the result, including an unsuccessful certification

output_dir = fullfile(project_root,'results');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

save(fullfile(output_dir,'theorem6.mat'), ...
    'p','Nx','Nt','y_0','y_f_base','args','source_witness','source_Nx', ...
    'certified','J','support','pairing','err','discretization_bound', ...
    'preparation_seconds','certification_seconds', 'bounds');

disp('Certified non-reachability:');
disp(certified);
disp('Final certification enclosure:');
disp(J);
disp('Certification time in seconds:');
disp(certification_seconds);
