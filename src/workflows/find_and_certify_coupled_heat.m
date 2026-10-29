function [J_opt,p_opt,u_opt] = find_and_certify_coupled_heat(y_0,y_0_base,y_f,y_f_base,p_0,u_0,args,args_float,target_kind, target_radius)
% Search for a coupled-heat witness and evaluate its verified certificate.
% Nonzero initial states are supported.
% Nodal vectors concatenate the first component and then the second.
% y_0 contains both nodal components.
% y_f contains the first nodal component, or both components.
% y_0_base and y_f_base describe both continuous components for verification.
% args and args_float are returned by init_coupled_heat_fem.
% target_kind: 1 = partial target (default), 2 = L2 ball
%
% J_opt is the verified functional enclosure at the returned witness p_opt.
% u_opt is the final primal iterate, or u_0 if the initial witness suffices.
% No optimality guarantee is provided.
% The returned witness uses nodal normalization, not exact spline L2 normalization.


% y_f_base describes both components, as required by the verifier.
if nargin < 9
    target_kind = 1;
end
validateattributes(target_kind, {'double'}, ...
    {'scalar','integer','>=',1,'<=',2});

% The radius is required for a ball and unused for a partial target.
if target_kind == 2
    assert(nargin >= 10, 'Provide target_radius for an L2 ball.');
    assert(isa(target_radius,'intval') && isequal(size(target_radius),[1,1]) ...
        && isreal(target_radius) && isfinite(inf(target_radius)) ...
        && isfinite(sup(target_radius)) && inf(target_radius) >= 0, ...
        'target_radius must be a finite, nonnegative scalar intval.');
else
    target_radius = intval(0);
end

Nt = args{2};
norm_L_T = args{7};
euler_powers = args_float{8};

if any(~isfinite(euler_powers), "all")
    error("This workflow requires a finite floating propagation table. " ...
        + "Initialize with init_coupled_heat_fem(Nx,Nt,intval_T,true).");
end

Nx = args_float{11};
Dx = 1/Nx;

coeffs = args{13};
a = coeffs(1,1); b = coeffs(1,2); c = coeffs(2,1); d = coeffs(2,2);
k1 = args{12}(1);
k2 = args{12}(2);

lambda_max_S_1 = (a+d+sqrt((a-d)^2 + (b+c)^2))/intval('2');
lambda_max_S_2 = a*k1 + d*k2 + sqrt((a*k1-d*k2)^2 + (b*k1 + c*k2)^2);
lambda_max_S_3 = a*k1 + d*k2 + sqrt((a*k1-d*k2)^2 + (c*k1 + b*k2)^2);


if ~(lambda_max_S_1 < intval('pi')^2*min(k1,k2) && max(lambda_max_S_2,lambda_max_S_3)  < intval('pi')^2*min(k1^2,k2^2))
    error("The coefficients of A do not satisfy the constraints to have bounded discretisation errors.")
end

if length(y_f) == Nx-1
    y_f_complete = zeros(2*Nx-2,1); y_f_complete(1:Nx-1) = y_f;
elseif length(y_f) == 2*Nx-2
    y_f_complete = y_f(:);
else
    error("Data y_f1 not understood.")
end

fl_norm_y_0 = norm(y_0)/sqrt(Nx);
fl_norm_y_f = norm(y_f)/sqrt(Nx);

y_tilde_T_float = y_f_complete - euler_powers(:,:,Nt+1)*y_0;

%% Certify the initial witness

p_init = y_tilde_T_float;
if target_kind == 1
    p_init(Nx:2*Nx-2) = 0;
end


% Normalize using the nodal norm.
% This does not impose an exact unit L2 norm on the continuous spline.

initial_norm = norm(p_init);
assert(isfinite(initial_norm), ...
    'The initial witness contains non-finite values.');

if initial_norm > 0
    p_init = p_init/(initial_norm*sqrt(Dx));
end



[target_not_reachable,J_init] = certify_coupled_heat_nonreachability(p_init,y_0_base,y_f_base,args,target_kind, target_radius);



disp("J(p_init) = ");
disp(J_init)
disp("Non-reachable : " + string(target_not_reachable))

if ~target_not_reachable
%% Floating witness search

if norm_L_T <= 0
    % Heuristic step-size scale when the control-map norm is unavailable.
    norm_L_T = 1;
end

gain = 1000;
sigma = gain/norm_L_T;
tau = 1/(gain*norm_L_T);
theta = 1;

if target_kind == 1
    prox_F_star = @(p,s)fl_prox_partial_target(p,s,y_tilde_T_float);
    target_radius_fl = NaN; % Unused for the partial target.

elseif target_kind == 2
    target_radius_fl = mid(target_radius); % Floating search approximation.
    fl_norm_y_f = fl_norm_y_f + target_radius_fl;
    prox_F_star = @(p,s)fl_prox_ball_target(p,s,y_tilde_T_float,target_radius_fl);

else
    error("Unknown target kind.");
end



n_iter = 250;
step_tolerance = 1e-5;
continue_search = @(x,y,x_old,y_old,i) should_continue_search(x,y,x_old,y_old, i, n_iter, step_tolerance,fl_norm_y_0, fl_norm_y_f,y_tilde_T_float,args_float, target_kind, target_radius_fl);



prox_G = @(u,t)fl_project_box_controls(u,tau,args_float);

K = @(u)fl_apply_control_map(u,args_float);
K_t = @(p)fl_apply_control_adjoint(p,args_float);


minus_p_0 = -p_0;

[u_opt,minus_p_opt] = chambolle_pock(prox_F_star,prox_G,K,K_t,sigma,tau,theta,u_0,minus_p_0,continue_search);


% Normalize using the nodal norm.
% This does not impose an exact unit L2 norm on the continuous spline.

witness_norm = norm(minus_p_opt);
if ~isfinite(witness_norm) || witness_norm == 0
    error('CAPSearch:InvalidWitness', ...
        'The search returned a zero or non-finite witness.');
end

% Switch back from the Chambolle-Pock sign convention and normalize.
p_opt = -minus_p_opt/(witness_norm*sqrt(Dx));


%% Certify the returned witness


[target_not_reachable,J_opt] = certify_coupled_heat_nonreachability( ...
    p_opt,y_0_base,y_f_base,args,target_kind, target_radius);


disp("J(p_opt) = ");
disp(J_opt);
disp("Certified non-reachability: " + string(target_not_reachable));
else
J_opt = J_init;
p_opt = p_init;
u_opt = u_0;
end


end


%% Local functions

function bool = should_continue_search(~,y,~,y_old, i, n_iter, step_tolerance,norm_y_0,norm_y_f,y_tilde_T_float,args_float, target_kind, target_radius_fl)
% Floating-search stopping heuristics.
% Return true to continue and false to stop.
% Final certification is performed separately.

Nx = args_float{11};

current_norm = norm(y);
if current_norm == 0
    y_unit = y;
else
    y_unit = y/current_norm;
end
y_norm = y_unit*sqrt(Nx);

previous_norm = norm(y_old);
if previous_norm == 0
    y_old_unit = y_old;
else
    y_old_unit = y_old/previous_norm;
end

normalized_step = norm(y_unit-y_old_unit);

% Floating diagnostics only; certification is performed after the search.
if i == 0
    disp("% ; J_diag_float ; J_diag_float + estimated error ; normalized step / tolerance")
end


% This stopping diagnostic integrates the support of the masked P1 field,
% while retaining nodal target terms.
% It need not equal the nodal functional associated with the search.
% Non-reachability is decided only by the subsequent verified evaluation.
control_field = fl_apply_control_adjoint(-y_norm,args_float);

if target_kind == 1
    J_disc = fl_p1_box_field_support(control_field,args_float) ...
        + fl_partial_target_term(y_norm,y_tilde_T_float);

elseif target_kind == 2
    J_disc = fl_p1_box_field_support(control_field,args_float) ...
        + fl_ball_target_term(y_norm,y_tilde_T_float,target_radius_fl);

else
    error("Unknown target kind.");
end

% Print periodic floating diagnostics.
if rem(20*i,n_iter) == 0
    error_estimate = estimate_discretization_error(y_norm,norm_y_0,norm_y_f,args_float);

    fprintf('%3.0f%% ; %.6e ; %.6e ; %.3e\n', ...
        100*i/n_iter,J_disc,J_disc+error_estimate, ...
        normalized_step/step_tolerance);
end

if i <= 2
    bool = 1;
else
    if i > n_iter 
        bool = 0;
        disp("100%");
        disp("Max number of iterations reached : " + num2str(i));
    else
        if J_disc < -1.2*estimate_discretization_error(y_norm,norm_y_0,norm_y_f,args_float)
            bool = 0;
            disp("100%");
            disp("Number of iterations : " + num2str(i));
            disp("Reason for stopping: floating functional below the estimated error margin.");
        else
            if normalized_step <= step_tolerance
                bool = 0;
                disp("100%");
                disp("Number of iterations : " + num2str(i));
                disp("Reason for stopping: normalized iterates stopped moving.");
            else
                bool = 1;
            end
        end
    end
end

end


function err = estimate_discretization_error(y_norm,norm_y_0,norm_y_f,args_float)
% Estimate discretization errors for floating-search stopping heuristics.
% y_norm contains the two scaled nodal witness components.
% norm_y_0 and norm_y_f are floating norm estimates.
% For a ball target, norm_y_f also includes the radius.
% This function does not provide a verified error bound.
% Final certification uses iv_coupled_heat_discretization_bound instead.

T = args_float{1};
Dt = args_float{3};
A = args_float{4};
coeffs_coerc_cont = args_float{6};
h = args_float{10};
Nx = args_float{11};
kappa_1 = args_float{12}(1);
kappa_2 = args_float{12}(2);

reaction_matrix = args_float{13};
a = reaction_matrix(1,1); b = reaction_matrix(1,2);
c = reaction_matrix(2,1); d = reaction_matrix(2,2);

M = max(abs(args_float{14})); % Amplitude bound from the control constraints.

% Floating P1 mass matrix and negative stiffness matrix.
int_range = linspace(1,Nx-1,Nx-1);
mat_D_diag =( 1 + 0.5*cos(pi/Nx*int_range))/((3/2)*Nx);

mat_V = sqrt(2/Nx)*sin(pi/Nx*int_range.*transpose(int_range));

mat_M = mat_V*diag(mat_D_diag)*mat_V;

toeplitz_K =  zeros(Nx-1,1);
toeplitz_K(1) = -2*Nx;
toeplitz_K(2) = 1*Nx;

mat_K = toeplitz(toeplitz_K);

% Heuristic derivative estimates from the discrete spatial operator.
norm_d2_p_f_1 = ...
    sqrt(h)*norm(mat_M\mat_K*y_norm(1:Nx-1));

norm_d2_p_f_2 = ...
    sqrt(h)*norm(mat_M\mat_K*y_norm(Nx:2*Nx-2));

% Lumped-mass estimate of the discrete adjoint-generator norm.
norm_A_star_p_f = sqrt(h)*norm(A'*y_norm);


% Bound on admissible forcings for pointwise bounded controls.
control_intervals = args_float{15};
omega_measure = sum(sup(control_intervals) - inf(control_intervals));
M_BU = M*sqrt(omega_measure);

lambda_max_S_2 = a*kappa_1 + d*kappa_2 ...
    + sqrt((a*kappa_1-d*kappa_2)^2 + (b*kappa_1+c*kappa_2)^2);

lambda_max_S_3 = a*kappa_1 + d*kappa_2 ...
    + sqrt((a*kappa_1-d*kappa_2)^2 + (c*kappa_1+b*kappa_2)^2);

delta_2 = min(kappa_1^2,kappa_2^2) - max(0,lambda_max_S_2)/pi^2;
delta_3 = min(kappa_1^2,kappa_2^2) - max(0,lambda_max_S_3)/pi^2;

C_0 = (1/sqrt(delta_2) + 1/sqrt(delta_3))/pi;

C_alpha = 2 + 2/sqrt(3);
a_0 = coeffs_coerc_cont(1);
a_1 = coeffs_coerc_cont(2);

% C_1 depends on whether A is self-adjoint or not
if b == c
    C_1 = a_1^(3/2)*C_0^2/(4*sqrt(a_0));
else
    C_1 = a_1^2*C_0^2/a_0;
end

C_2 = C_1*(7+4*log(2)*a_1/a_0+C_alpha);


C_3 = C_alpha*a_1/(2*exp(1)*a_0);

err_1 = (0.5*Dt*M_BU*T ...
    + (norm_y_0 + M_BU*T)*(C_2*h^2 + C_3*Dt)) ...
    * norm_A_star_p_f;

interpolation_error_bound = ...
    sqrt(norm_d2_p_f_1^2 ...
    + norm_d2_p_f_2^2)*h^2/pi^2;

err_2 = ((norm_y_0 + M_BU*T)*C_alpha + norm_y_f) ...
    * interpolation_error_bound;

err = err_1 + err_2;
end
