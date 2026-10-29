% Proposition 20: unregularized witness search for illustration only.
% This script does not call a non-reachability verifier or produce a proof.
% It reuses the saved coarse model instead of the current initializer settings.

project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(project_root,'src')));
witness_dir = fullfile(project_root,'data','witnesses');

%% Reuse the coarse tent configuration

regularized = load(fullfile(witness_dir,'proposition20_coarse.mat'));
args = regularized.args;
Nx = args{11};
Nt = args{2};
y_f_base = regularized.y_f_base;
y_f = mid(y_f_base);

assert(Nx <= 100 && Nt <= 2000, ...
    'Use the saved coarse search configuration, not a fine certification mesh.');
assert(mod(Nx,5) == 0 && length(y_f) == Nx-1, ...
    'The exact P1 tent must have its vertex on the mesh.');
assert(inf(args{1}) == 1 && sup(args{1}) == 1, ...
    'Proposition 20 requires final time T = 1.');

%% Recover floating arguments and check both control intervals

args_float = args;
for k = 1:length(args_float)
    if isa(args_float{k},'intval')
        args_float{k} = mid(args_float{k});
    end
end
args_float{13} = args{13};

x = (1:Nx-1)'/Nx;
control_mask = false(Nx-1,1);
for k = 1:length(args{13})
    control_mask = control_mask | ...
        (x >= inf(args{13}(k)) & x <= sup(args{13}(k)));
end
assert(isequal(args_float{5},diag(double(control_mask))), ...
    'The saved control mask must include every control interval.');

args_float{8} = fl_euler_power_table( ...
    Nt,args_float{3},args_float{4});

%% Search parameters: lambda = 0, with the same initial guesses

lambda_reg = 0;
gain = 10000;
max_iterations = 1000;
direction_tolerance = 1e-10;

norm_L_T = args_float{7};
if norm_L_T <= 0
    norm_L_T = 1;
end
sigma = gain/norm_L_T;
tau = 1/(gain*norm_L_T);
theta = 1;

p_0 = y_f/(norm(y_f)*sqrt(1/Nx));
u_0 = zeros(Nx-1,Nt);

prox_F_star = @(q,s) fl_prox_heat_regularized_target( ...
    q,s,y_f,lambda_reg,args_float);
prox_G = @(v,t) fl_project_positive_controls(v,t,args_float);
K = @(v) fl_apply_control_map(v,args_float);
K_star = @(q) fl_apply_control_adjoint(q,args_float);

keep_going = @(v,q,v_old,q_old,k) continue_search( ...
    q,q_old,k,max_iterations,direction_tolerance,args_float,y_f);

%% Chambolle-Pock only: no certification is performed

started = tic;
[u,minus_p] = chambolle_pock( ...
    prox_F_star,prox_G,K,K_star,sigma,tau,theta,u_0,-p_0,keep_going);
search_seconds = toc(started);

assert(norm(minus_p) > 0,'The search returned a zero witness.');
p = -minus_p/(norm(minus_p)*sqrt(1/Nx));
J_search = fl_heat_positive_control_support(p,args_float) - p'*y_f/Nx;

%% Save separately from the regularized witness and certified results

certification_performed = false;
proof_status = 'not_certified';
output_file = fullfile(witness_dir,'proposition20_unregularized.mat');

save(output_file,'p','u','Nx','Nt','args','y_f_base', ...
    'lambda_reg','gain','max_iterations','direction_tolerance', ...
    'J_search','search_seconds','certification_performed','proof_status');

disp('Unregularized floating functional (not a certified bound):');
disp(J_search);
disp('Unregularized witness saved to:');
disp(output_file);

function keep_going = continue_search(q,q_old,k,max_iterations,tolerance,args_float,target)
% Match the direction-based stopping rule used by the regularized search.
if norm(q) > 0
    p = -q/(norm(q)*sqrt(args_float{10}));
else
    p = zeros(size(q));
end
if mod(k,50) == 0
    value = fl_heat_positive_control_support(p,args_float) ...
        - p'*target*args_float{10};
    fprintf('Iteration %d: floating J = %+.8e\n',k,value);
end
if k >= max_iterations
    fprintf('Stopping: maximum iteration count reached.\n');
    keep_going = false;
elseif k >= 3 && norm(q) > 0 && norm(q_old) > 0 ...
        && norm(q/norm(q)-q_old/norm(q_old)) <= tolerance
    fprintf('Stopping: normalized witness direction stabilized.\n');
    keep_going = false;
else
    keep_going = true;
end
end
