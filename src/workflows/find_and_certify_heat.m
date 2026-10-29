function [J_opt,p_opt_reg,u_opt] = find_and_certify_heat(y_f,y_f_base,p_0,u_0,args,args_float,lambda_reg)
% Search for a witness and evaluate its non-reachability certificate.
% This scalar workflow supports only the zero initial state.
%
% Inputs:
% - y_f: floating target values at the interior mesh nodes.
% - y_f_base: continuous target coefficients, interpreted by their length:
%   Nx-1 for P1; 2*Nx for cubic Hermite; otherwise sine-series coefficients.
%   Cubic Hermite data contain interior values followed by Dx times
%   the derivatives at all nodes, including the endpoints.
% - p_0, u_0: initial dual and primal iterates for the floating search.
% - args, args_float: data returned by init_heat_fem.
% - lambda_reg: floating-search regularization parameter; defaults to zero.
%
% Outputs:
% - J_opt: verified functional enclosure at the returned witness.
% - p_opt_reg: returned witness, scaled so sqrt(Dx)*norm(p_opt_reg) = 1.
% - u_opt: final primal iterate, or u_0 if the initial witness suffices.
% No optimality guarantee is provided.



%% Extraction of useful parameters


norm_L_T = args{7};
Nx = args_float{11};
Dx = 1/Nx;
euler_powers = args_float{8};

if nargin < 7
    lambda_reg = 0;
end

% This scalar workflow supports only the zero initial state.
y_0 = zeros(Nx-1,1);
y_tilde_T_float = y_f;


%% Certify the normalized target-shaped witness

initial_norm = norm(y_f);
assert(isfinite(initial_norm) && initial_norm > 0, ...
    'Cannot normalize the target-shaped initial witness.');

p_init = y_f/(initial_norm*sqrt(Dx));

[target_not_reachable,J_init] = ...
    certify_heat_nonreachability(p_init,y_0,y_f_base,args);

disp("J(p_init) = ");
disp(J_init)



if ~target_not_reachable %If p_f = y_f didn't show negativity

if any(~isfinite(euler_powers), "all")
    error("Witness search requires a finite floating propagation table. " ...
        + "Initialize with init_heat_fem(Nx,Nt,intval_T,true).");
end

%% Non-guaranteed search for a negative point (withouth Intlab)


% Values affecting the convergence (rate) of the Chambolle-Pock algorithm
gain = 10000; % Balance the primal and dual step sizes.

if norm_L_T < 0 % The initializer did not compute the control-map norm.
    norm_L_T = 1;
end

sigma = gain/norm_L_T;
tau = 1/(gain*norm_L_T);
theta = 1;

n_iter = 1000; %Max number of iterations. ~100 is a good number
step_tolerance = 1e-10; % Stopping condition if the renormalised minimiser doesn't more than step_tolerance

%Summary of the stopping conditions:
continue_search = @(x,y,x_old,y_old,i) should_continue_search(x,y,x_old,y_old, i, n_iter, step_tolerance,y_tilde_T_float,y_0,y_f,args_float);


%Proximal operators:
prox_F_star = @(p,s)fl_prox_heat_regularized_target(p,s,y_tilde_T_float, lambda_reg,args_float);
prox_G = @(u,t)fl_project_positive_controls(u,tau,args_float); %This is the projection on the constraints, here 0 <= u <= M

%Control-to-state operator and its dual
K = @(u)fl_apply_control_map(u,args_float);
K_t = @(p)fl_apply_control_adjoint(p,args_float);


minus_p_0 = -p_0; % Switching from control conventions to CP2011 conventions

%Chambolle-Pock algorithm :
[u_opt,minus_p_opt] = chambolle_pock(prox_F_star,prox_G,K,K_t,sigma,tau,theta,u_0,minus_p_0,continue_search);

witness_norm = norm(minus_p_opt);
assert(isfinite(witness_norm) && witness_norm > 0, ...
    'The floating search returned a zero or non-finite witness.');

% Normalize using the nodal norm.
% This does not impose an exact unit L2 norm on the continuous spline.
p_opt = -minus_p_opt/(witness_norm*sqrt(Dx));


%% Guaranteed negativity check (with Intlab)

% The regularization is included in the floating search.
% Certify the returned witness without additional smoothing.
p_opt_reg = p_opt;

%Computation of verified value of J
[~,J_opt] = certify_heat_nonreachability(p_opt_reg,y_0,y_f_base,args);

disp("J(p_opt) = ");
disp(J_opt);


else % If the case y_f = p_f was enough to show negativity
    J_opt = J_init;
    p_opt_reg = p_init;
    u_opt = u_0;
end

end




%% Annex functions



function bool = should_continue_search(~,y_i,~,y_im1, i, n_iter, step_tolerance,y_tilde_T_float,y_0,y_f,args_float)
% Stopping condition for the Chambolle-Pock witness search.
% Takes x, x_old, y, y_old, i, n_iter, step_tolerance, y_tilde_T_float, y_0, y_f, args_float :
% - i is the iteration number
% - n_iter is the maximal interation number
% - x_im1, x_i are the primal values at iterations i-1, i. Useless here
% - y_im1,y_i are the dual values at iterations i-1, i.
% - y_tilde_T_float = y_f - S_Ty_0 the numerical target
% - y_0 the numerical initial state
% - y_f the numerical final state
% - args_float the values returned by init_heat_fem
%
% After the first three iterations, stop when:
% - the iteration limit is exceeded;
% - the floating functional is below -1.2 times the estimated error;
% - the Euclidean distance between the scaled iterates is at most step_tolerance.
% These are search heuristics, not certification conditions.
%
% Returns 1 if the algorithm can continue, 0 if it must stop
% Prints interesting values during the iterations

%% Parameters needed
Nx = args_float{11};

%Renormalisation of y_i
if norm(y_i) == 0 
    y_norm_i = y_i;
else
    y_norm_i = y_i/norm(y_i) * sqrt(Nx);
end

%Renormalisation of y_im1
if norm(y_im1) == 0 
    y_norm_im1 = y_im1;
else
    y_norm_im1 = y_im1/norm(y_im1) * sqrt(Nx);
end

% Floating diagnostics only; certification is performed after the search.
if i == 0
    disp("% ; J_float ; J_float + estimated error ; normalized step / tolerance")
end

% Print every 5% of the iteration limit.
if rem(20*i,n_iter) == 0
    % Convert from the Chambolle-Pock sign convention.
    J_float = fl_heat_positive_control_support(-y_norm_i,args_float) ...
        + transpose(y_norm_i)*y_tilde_T_float/Nx;

    error_estimate = estimate_discretization_error(y_norm_i,y_0,y_f,args_float);
    step_ratio = norm(y_norm_i-y_norm_im1)/step_tolerance;

    fprintf('%3.0f%% ; %.6e ; %.6e ; %.3e\n', ...
        100*i/n_iter,J_float,J_float+error_estimate,step_ratio);
end

%% Stopping condition

if i <= 2
    bool = 1;
elseif i > n_iter 
    bool = 0;
    disp("Max number of iterations reached : " + num2str(i));
elseif (fl_heat_positive_control_support(-y_norm_i, args_float) + transpose(y_norm_i)*y_tilde_T_float/Nx) < -1.2*estimate_discretization_error(y_norm_i,y_0,y_f,args_float)
    bool = 0;
    disp("Reason of stopping : J_d(p_opt) sufficiently negative")
elseif norm(y_norm_i - y_norm_im1) <= step_tolerance 
    bool = 0;
    disp("Reason of stopping : the algorithm stopped moving.")
else
    bool = 1;
end
end


function err = estimate_discretization_error(y_norm,y_0,y_f,args_float)
% Estimate discretization errors for floating-search stopping heuristics.
% y_norm is a scaled nodal witness; y_0 and y_f are nodal state vectors.
% args_float contains the floating data returned by init_heat_fem.
% This function does not provide a verified error bound.
% Final certification uses iv_heat_discretization_bound instead.

% Parameters
T = args_float{1};
Dt = args_float{3};
alpha = args_float{6};
M = args_float{9};
Dx = args_float{10};
Nx = args_float{11};
kappa = args_float{12};

% Heuristic derivative and interpolation estimates from the discrete operator.
norm_A_star_p_f = norm(args_float{4}*y_norm)/sqrt(Nx);
norm_d2_p_f = norm_A_star_p_f/kappa;
interpolation_error_bound = Dx^2/(pi^2)*norm_d2_p_f;

% y_0 always assumed to be P1
norm_y_0 = norm(y_0)*sqrt(Dx);

% Lumped-mass approximation of the target L2 norm from its nodal values.
norm_y_f = norm(y_f)/sqrt(Nx);

% Bound on admissible forcings for positive pointwise controls.
control_intervals = args_float{13};
omega_measure = sum(sup(control_intervals) - inf(control_intervals));
M_BU = M*sqrt(omega_measure);

% Floating constants used by the heuristic error estimate.
C_0 = 2/(pi*kappa);
C_alpha = 2 + 2/sqrt(3);
a_0 = kappa;
a_1 = kappa;
C_1 = a_1^(3/2)*C_0^2/(4*sqrt(a_0));
C_2 = C_1*(7+4*log(2)/cos(alpha)+C_alpha);
C_3 = C_alpha/(2*exp(1)*cos(alpha));

%Error computations 
err_1 = (0.5*M_BU*T*Dt ...
    + (norm_y_0 + M_BU*T)*(C_2*Dx^2 + C_3*Dt)) ...
    * norm_A_star_p_f;

err_2 = ((norm_y_0 + M_BU*T)*C_alpha + norm_y_f) ...
    * interpolation_error_bound;


err = err_1 + err_2;
end
