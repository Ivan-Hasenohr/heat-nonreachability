function [err,bounds] = iv_heat_discretization_bound(p_f_base,y_0,norm_y_f,args,type_constraints,verbose)
% Bound the discretization error for the scalar heat equation.
%
% Inputs:
% - p_f_base: cubic Hermite coefficient vector of length 2*Nx,
%   with zero endpoint values.
% - y_0: initial state; only the zero initial state is supported here.
% - norm_y_f: verified target-norm bound for the target-pairing error.
%   Pass zero when that pairing is evaluated directly on the spline witness.
% - args: verified data from init_heat_fem.
% - type_constraints: 1 = pointwise bounds; 2 = L2 bounds (default).
% - verbose: optional diagnostic output; false by default.
%
% Outputs:
% - err: upper bound on the total discretization error.
% - bounds: time, space and total discretization bounds.
% Interval uncertainty in the discrete functional is recorded by the caller.

% Optional arguments: control constraints and diagnostic output.
if nargin < 5 || isempty(type_constraints)
    type_constraints = 2;
end
if nargin < 6 || isempty(verbose)
    verbose = false;
end

validateattributes(type_constraints, {'double'}, ...
    {'scalar','integer','>=',1,'<=',2});

%Parameters
T = args{1};
Dt = args{3};

alpha = args{6};
M = args{9}; 
h = args{10};
kappa = args{12};
control_intervals = args{13};


% ||p_f''||
norm_d2_p_f = iv_cubic_hermite_d2_norm(p_f_base);

% ||A^* p_opt||
norm_A_star_p_f = kappa*norm_d2_p_f;

% Bound on ||p_f - p_fh||
interpolation_error_bound = h^2/intval('pi')^2 * norm_d2_p_f;

% Only the zero initial state is supported by this scalar implementation.
initial_data = intval(y_0);
assert(all(inf(initial_data(:)) == 0) ...
    && all(sup(initial_data(:)) == 0), ...
    'Only the zero initial state is supported.');

norm_y_0 = intval(0);

% Bound on admissible forcings, distinct from the operator norm of B.
if type_constraints == 1
    omega_measure = intval(0);

    for i = 1:length(control_intervals)
        omega_measure = omega_measure ...
            + intval(sup(control_intervals(i))) ...
            - intval(inf(control_intervals(i)));
    end

    M_BU = M*sqrt(omega_measure);
else
    % Restriction to the control domain has L2 operator norm at most one.
    M_BU = M;
end

%% Constants

C_0 = 2/(intval('pi')*kappa);
C_alpha = 2+2/intval('sqrt3');
a_0 = kappa;
a_1 = kappa;

% The scalar Dirichlet diffusion operator is self-adjoint.
C_1 = a_1^(3/2)*C_0^2/(4*sqrt(a_0));

C_2 = C_1*(7+4*intval('log2')/cos(alpha)+C_alpha);
C_3 = C_alpha/(2*exp(intval(1))*cos(alpha));

%Computation of errors
err_1 = (0.5*M_BU*T*Dt + (norm_y_0 + M_BU*T)*(C_2*h^2+C_3*Dt))*norm_A_star_p_f;
err_2 = ((norm_y_0 + M_BU*T)*C_alpha + norm_y_f)*interpolation_error_bound;

if verbose
    disp("Time error bound: " + string(sup( ...
        Dt*norm_A_star_p_f*(0.5*M_BU*T ...
        + C_3*(norm_y_0 + M_BU*T)))));

    disp("Space error bound: " + string(sup( ...
        h^2*norm_A_star_p_f*C_2*(norm_y_0 + M_BU*T) ...
        + err_2)));
end


err = sup(err_1 + err_2);

% Optional diagnostic bounds; the certified total above is unchanged.
if nargout > 1
    bounds.time = sup( ...
        Dt*norm_A_star_p_f*(0.5*M_BU*T ...
        + C_3*(norm_y_0 + M_BU*T)));

    bounds.space = sup( ...
        h^2*norm_A_star_p_f*C_2*(norm_y_0 + M_BU*T) ...
        + err_2);

    bounds.discretization = err;
end