function [err,bounds] = iv_coupled_heat_discretization_bound(p_opt_base,norm_y_0,norm_y_f,args,verbose)
% Bound the discretization error for the coupled heat equation.
%
% Inputs:
% - p_opt_base: two concatenated cubic Hermite coefficient vectors,
%   each of length 2*Nx, with zero endpoint values.
% - norm_y_0: verified upper bound on the product L2 norm of y_0.
% - norm_y_f: verified target-norm bound for the target-pairing error.
%   Pass zero when that pairing is evaluated directly on the spline witness.
% - args: verified data from init_coupled_heat_fem.
% - verbose: optional diagnostic output; false by default.
%
% Outputs:
% - err: upper bound on the total discretization error.
% - bounds: time, space and total discretization bounds.
% Interval uncertainty in the discrete functional is recorded by the caller.

% Disable diagnostic output by default.
if nargin < 5 || isempty(verbose)
    verbose = false;
end

assert(isa(args{1},'intval') ...
    && isa(args{3},'intval') ...
    && isa(args{10},'intval'), ...
    'T, Dt and h must be provided as verified intervals.');

T = args{1};
Dt = args{3};
form_bounds = args{6};
h = args{10};
Nx = args{11};
kappa_1 = args{12}(1);
kappa_2 = args{12}(2);

reaction_matrix = args{13};
a = reaction_matrix(1,1); b = reaction_matrix(1,2); c = reaction_matrix(2,1); d = reaction_matrix(2,2);

% Derive the amplitude bound from the actual control constraints.
M = max(abs(intval(args{14})));


% Validate the mesh contract before evaluating the certificate.
Nt_check = args{2};
validateattributes(Nt_check,{'double'}, ...
    {'real','finite','scalar','integer','>=',2});
validateattributes(Nx,{'double'}, ...
    {'real','finite','scalar','integer','>=',3});

for entry = [1,3,10]
    value = args{entry};
    assert(isequal(size(value),[1,1]) && isreal(value) ...
        && isfinite(inf(value)) && isfinite(sup(value)) && inf(value) > 0, ...
        'T, Dt and h must be finite positive scalar intervals.');
end

expected_Dt = args{1}/intval(Nt_check);
expected_h = intval(1)/intval(Nx);

assert(inf(args{3}) <= inf(expected_Dt) ...
    && sup(args{3}) >= sup(expected_Dt), ...
    'Dt must enclose T/Nt.');

assert(inf(args{10}) <= inf(expected_h) ...
    && sup(args{10}) >= sup(expected_h), ...
    'h must enclose 1/Nx.');

%%  Computation of norm of 2nd derivatives:

norm_d2_p_f_1 = iv_cubic_hermite_d2_norm(p_opt_base(1:2*Nx));
norm_d2_p_f_2 = iv_cubic_hermite_d2_norm(p_opt_base(2*Nx+1:4*Nx));

%% L2 Gram matrix in the cubic Hermite basis

intval_M = intval(zeros(2*Nx,2*Nx));
intval_M(1:Nx-1,1:Nx-1) = 78*eye(Nx-1) + intval('27')/intval('2')*(intval(diag(ones(Nx-2,1),1)) + diag(ones(Nx-2,1),-1));
intval_M(Nx:2*Nx,Nx:2*Nx) = 2*eye(Nx+1) -intval('3')/intval('4')*( diag(ones(Nx,1),1) + diag(ones(Nx,1),-1));
intval_M(Nx,Nx) = 1;intval_M(2*Nx,2*Nx) = 1;



intval_M(1:Nx-1,Nx:2*Nx-2) = intval('13')/intval('4')*diag(ones(Nx-1,1)) ;intval_M(1:Nx-1,Nx+2:2*Nx) = intval_M(1:Nx-1,Nx+2:2*Nx) - intval('13')/intval('4')*diag(ones(Nx-1,1));

intval_M(Nx:2*Nx,1:Nx-1) = transpose(intval_M(1:Nx-1,Nx:2*Nx));

intval_M = intval_M/(intval('105')*intval(Nx));

%% Negative H1 seminorm Gram matrix in the cubic Hermite basis

intval_K = intval(zeros(2*Nx,2*Nx));
intval_K(1:Nx-1,1:Nx-1) = 36*eye(Nx-1) -18*( diag(ones(Nx-2,1),1) + diag(ones(Nx-2,1),-1) );
intval_K(Nx:2*Nx,Nx:2*Nx) = 4*eye(Nx+1) -(1/intval('2'))*( diag(ones(Nx,1),1) + diag(ones(Nx,1),-1) );
intval_K(Nx,Nx) = 2;intval_K(2*Nx,2*Nx) = 2;

intval_K(1:Nx-1,Nx:2*Nx-2) = -(intval('3')/intval('2'))* diag(ones(Nx-1,1));intval_K(1:Nx-1,Nx+2:2*Nx) = intval_K(1:Nx-1,Nx+2:2*Nx) + (intval('3')/intval('2'))* diag(ones(Nx-1,1));
intval_K(Nx:2*Nx,1:Nx-1) = transpose(intval_K(1:Nx-1,Nx:2*Nx));


intval_K = -intval_K*intval(Nx)/intval('15');

%% Squared L2 norms of the first derivatives

norm_d1_p_f_1_squared = -transpose(p_opt_base(1:2*Nx))*intval_K*p_opt_base(1:2*Nx);
norm_d1_p_f_2_squared = -transpose(p_opt_base(2*Nx+1:4*Nx))*intval_K*p_opt_base(2*Nx+1:4*Nx);

%% Computation of ||A^*p_f||

% Negative derivative pairing: -<p_1',p_2'>:
negative_d1_pairing = transpose(p_opt_base(1:2*Nx))*intval_K*p_opt_base(2*Nx+1:4*Nx); 

% Reaction components of A^* p_f: the coupling matrix is transposed.
adjoint_reaction_1 = a*p_opt_base(1:2*Nx) + c*p_opt_base(2*Nx+1:4*Nx);
adjoint_reaction_2 = b*p_opt_base(1:2*Nx) + d*p_opt_base(2*Nx+1:4*Nx);

norm_adjoint_reaction_1_squared = transpose(adjoint_reaction_1)*intval_M*adjoint_reaction_1; 
norm_adjoint_reaction_2_squared = transpose(adjoint_reaction_2)*intval_M*adjoint_reaction_2;

norm_A_star_p_f = sqrt(kappa_1^2*norm_d2_p_f_1^2 + kappa_2^2*norm_d2_p_f_2^2 + norm_adjoint_reaction_1_squared + norm_adjoint_reaction_2_squared + 2*(c*kappa_1 + b*kappa_2)*negative_d1_pairing -2*a*kappa_1*norm_d1_p_f_1_squared - 2*d*kappa_2*norm_d1_p_f_2_squared);

% Bound on admissible forcings for pointwise bounded controls.
control_intervals = args{15};
omega_measure = intval(0);

for i = 1:length(control_intervals)
    omega_measure = omega_measure + intval(sup(control_intervals(i))) - intval(inf(control_intervals(i)));
end

M_BU = M*sqrt(omega_measure);

%% Constants 

% Regularity constant accounting for the coupling.
lambda_max_S_2 = a*kappa_1 + d*kappa_2 + sqrt((a*kappa_1-d*kappa_2)^2 + (b*kappa_1+c*kappa_2)^2);

lambda_max_S_3 = a*kappa_1 + d*kappa_2 + sqrt((a*kappa_1-d*kappa_2)^2 + (c*kappa_1+b*kappa_2)^2);

delta_2 = min(kappa_1^2,kappa_2^2) - max(intval(0),lambda_max_S_2)/intval('pi')^2;

delta_3 = min(kappa_1^2,kappa_2^2) - max(intval(0),lambda_max_S_3)/intval('pi')^2;

if ~(inf(delta_2) > 0 && inf(delta_3) > 0)
    error("The regularity denominators are not certified positive.")
end

C_0 = (1/sqrt(delta_2) + 1/sqrt(delta_3))/intval('pi');


C_alpha = 2 + 2/intval('sqrt3');
a_0 = form_bounds(1);
a_1 = form_bounds(2);



% Certify symmetry only for equal point coefficients.
is_self_adjoint = inf(b) == sup(b) ...
    && inf(c) == sup(c) ...
    && inf(b) == inf(c);

% C_1 depends on whether A is self-adjoint or not
if is_self_adjoint
    if verbose
        disp("Self-adjoint case");
    end
    C_1 = a_1^(3/2)*C_0^2/(4*sqrt(a_0));
else
    if verbose
        disp("Using the general constant: self-adjointness is not certified.");
    end
    C_1 = a_1^2*C_0^2/a_0;
end

C_2 = C_1*(7+4*intval('log2')*a_1/a_0+C_alpha);
C_3 = C_alpha*a_1/(2*exp(intval(1))*a_0);
%% Computation of discretisation errors

err_1 = (0.5*Dt*M_BU*T + (norm_y_0 + M_BU*T)*(C_2*h^2+C_3*Dt))*norm_A_star_p_f;

interpolation_error_bound = sqrt(norm_d2_p_f_1^2 + norm_d2_p_f_2^2) * h^2/intval('pi')^2;

err_2 = ((norm_y_0 + M_BU*T)*C_alpha + norm_y_f)*interpolation_error_bound;

if verbose
    disp("Time error bound: " + string(sup( ...
        Dt*norm_A_star_p_f*(0.5*M_BU*T ...
        + C_3*(norm_y_0 + M_BU*T)))));
    disp("Space error bound: " + string(sup( ...
        err_2 + h^2*norm_A_star_p_f*C_2*(norm_y_0 + M_BU*T))));
    disp("||A^*p_f|| = " + string(sup(norm_A_star_p_f)));
end
err = sup(err_1 + err_2);

% Optional diagnostic bounds; the certified total above is unchanged.
if nargout > 1
    bounds.time = sup( ...
        Dt*norm_A_star_p_f*(0.5*M_BU*T ...
        + C_3*(norm_y_0 + M_BU*T)));

    bounds.space = sup( ...
        err_2 + h^2*norm_A_star_p_f*C_2*(norm_y_0 + M_BU*T));

    bounds.discretization = err;
end

end
