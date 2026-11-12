function result = run_coupled_projected_certificate(example_name,Nx,Nt,run_name)
%RUN_COUPLED_PROJECTED_CERTIFICATE Certify a coupled heat example.
% example_name: 'ball', 'partial' or 'split'. Initialize INTLAB first.
% Nx and Nt select the mesh; run_name labels files under results/coupled_projected.
% The original FEM/Euler certifier constructs a natural cubic witness from p.
% Example: run_coupled_projected_certificate('split',3090,47000,'fine').

example_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(example_dir));
addpath(genpath(fullfile(project_root,'src')));

validateattributes(Nx,{'double'},{'scalar','real','finite','integer','>=',3});
validateattributes(Nt,{'double'},{'scalar','real','finite','integer','>=',2});
assert(ischar(run_name) && isrow(run_name) ...
    && ~isempty(regexp(run_name,'^[A-Za-z0-9_-]+$','once')), ...
    'run_name must be a nonempty filename label.');

%% Define the continuous problem and the sine witness

switch example_name
    case 'ball'
        T = intval(8)/25;
        control_intervals = infsup(0,0.5);
        control_bounds = intval([-1,2]);
        target_kind = 2;
        target_radius = intval(1)/5000;
        regions = [intval(0),intval(1)/2];
        target = iv_coupled_constant_control_endpoint(2,T,regions,intval(9)/5);
        coefficients = [.1,0; -1,.6; -1,-1; 0,.6];

    case 'partial'
        T = intval(1)/5;
        control_intervals = infsup(0,0.5);
        control_bounds = intval([-1,2]);
        target_kind = 1;
        target_radius = intval(0);
        target = [-intval(23)/10000,0; intval(1)/1000,0];
        coefficients = [-1,0; 6,0];

    case 'split'
        T = intval(3)/50;
        control_intervals = infsup(0,0.5);
        control_bounds = intval([-1,2]);
        target_kind = 2;
        target_radius = intval(0);
        quarter = intval(1)/4;
        half = intval(1)/2;
        % The two contributions are evaluated directly at c = 1/2.
        target = iv_coupled_constant_control_endpoint(4,T,[intval(0),quarter],intval(1)) ...
            + iv_coupled_constant_control_endpoint(4,T,[quarter,half],-half);
        coefficients = readmatrix(fullfile(project_root,'data','witnesses', ...
            'coupled_split_sine.csv'));

    otherwise
        error('Choose ball, partial or split.');
end

%% Initialize the verified FEM discretization

fprintf('\n%s: Nx = %d, Nt = %d (%s).\n',example_name,Nx,Nt,run_name);
started = tic;
args = init_coupled_heat_fem(Nx,Nt,T,false);

% The target formula and the auxiliary checks use this specific operator.
% Fail explicitly if the initializer's physical parameters have changed.
reaction = intval([-1,2;-2,-1]);
assert(isequal(inf(args{13}),inf(reaction)) ...
    && isequal(sup(args{13}),sup(reaction)) ...
    && isequal(inf(args{12}),[1;1.5]) && isequal(sup(args{12}),[1;1.5]), ...
    'The projected targets require the documented coupled operator.');

% All examples use the initializer's control region (0,1/2).
assert(isequal(inf(args{15}),0) && isequal(sup(args{15}),0.5), ...
    'These examples require the control region (0,1/2).');
args{14} = control_bounds;
args{9} = max(abs(control_bounds));
x = (1:Nx-1)'/Nx;

%% Sample the fixed sine witness on the chosen mesh

assert(Nx>size(coefficients,1), 'Use more cells than witness frequencies.');
values = sin(pi*x*(1:size(coefficients,1)))*coefficients;
p = [values(:,1);values(:,2)];
assert(norm(p)>0 && all(isfinite(p)), 'The sampled witness must be nonzero.');
% Joint nodal normalization; the continuous spline norm is checked separately.
p = p/(norm(p)*sqrt(1/Nx));
if target_kind == 1
    assert(all(p(Nx:end) == 0));
end

y_0_base = intval([0;0]);
y_f_base = [target(:,1);target(:,2)];
assert(length(y_f_base) ~= 2*(Nx-1) && length(y_f_base) ~= 4*Nx);

fprintf('Initialization finished after %.1f s. Starting certification.\n',toc(started));

%% Run the original certification function

[certified,J,bounds] = certify_coupled_heat_nonreachability( ...
    p,y_0_base,y_f_base,args,target_kind,target_radius);
seconds = toc(started);

%% Save the witness, error bounds and final enclosure

output_dir = fullfile(project_root,'results','coupled_projected',example_name);
if ~isfolder(output_dir)
    mkdir(output_dir);
end
output_file = fullfile(output_dir,run_name);
kappa = args{12};
form_bounds = args{6};
save([output_file '.mat'],'certified','J','bounds','p','coefficients', ...
    'Nx','Nt','T','y_0_base','y_f_base','target_kind','target_radius','seconds', ...
    'kappa','reaction','form_bounds','control_intervals','control_bounds');

result.example = example_name;
result.run_name = run_name;
result.Nx = Nx;
result.Nt = Nt;
result.certified = logical(certified);
result.J = [inf(J),sup(J)];
result.discrete_J = [inf(bounds.discrete_enclosure),sup(bounds.discrete_enclosure)];
result.time_error = bounds.time;
result.space_error = bounds.space;
result.roundoff = bounds.roundoff;
result.seconds = seconds;

fid = fopen([output_file '.json'],'w');
assert(fid >= 0);
fprintf(fid,'%s\n',jsonencode(result));
fclose(fid);

fprintf('Discrete J: [%.9g, %.9g].\n',result.discrete_J);
fprintf('Time error: %.9g; space error: %.9g; interval radius: %.9g.\n', ...
    result.time_error,result.space_error,result.roundoff);
fprintf('Final J: [%.9g, %.9g]; certified = %d; %.1f s.\n', ...
    result.J,result.certified,seconds);
end
