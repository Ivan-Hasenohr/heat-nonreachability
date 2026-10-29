function [target_not_reachable,J_opt,bounds] = certify_coupled_heat_nonreachability(p_opt,y_0_base,y_f_base,args,target_kind, target_radius)

% target_kind: 1 = partial target, 2 = L2 ball
if nargin < 5
    target_kind = 1;
end

validateattributes(target_kind, {'double'}, {'scalar','integer','>=',1,'<=',2});

% The radius is required for a ball and unused for a partial target.
if target_kind == 2
    assert(nargin >= 6, 'Provide target_radius for an L2 ball.');
    assert(isa(target_radius,'intval') && isequal(size(target_radius),[1,1]) ...
        && isreal(target_radius) && isfinite(inf(target_radius)) ...
        && isfinite(sup(target_radius)) && inf(target_radius) >= 0, ...
        'target_radius must be a finite, nonnegative scalar intval.');
else
    target_radius = intval(0);
end

Nx = args{11};

% Certification requires the verified initializer output.
assert(isa(args{1},'intval') ...
    && isa(args{3},'intval') ...
    && isa(args{4},'intval') ...
    && isa(args{10},'intval'), ...
    'Use verified args from init_coupled_heat_fem, not args_float.');


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


assert(isequal(size(args{4}),[2*Nx-2,2*Nx-2]), ...
    'The verified evolution matrix has an inconsistent size.');


% The fixed witness contains two concatenated nodal components.
validateattributes(p_opt,{'double'}, ...
    {'real','finite','column','numel',2*(Nx-1)});

% Both components use the same representation and coefficient count.
assert(size(y_0_base,2) == 1 && length(y_0_base) >= 2 ...
    && mod(length(y_0_base),2) == 0, ...
    'The initial state must contain two equally sized column components.');

assert(size(y_f_base,2) == 1 && length(y_f_base) >= 2 ...
    && mod(length(y_f_base),2) == 0, ...
    'The target data must contain two equally sized column components.');

p_opt_base_1 = iv_natural_cubic_from_nodes(p_opt(1:Nx-1));
p_opt_base_2 = iv_natural_cubic_from_nodes(p_opt(Nx:2*Nx-2));
p_opt_base = [p_opt_base_1 ; p_opt_base_2];


if target_kind == 1
    disp("Target of the form  Y_f = {(y_1,y) for y \in L^2(0,1)}")
elseif target_kind == 2
    disp("Target: L2 ball centered at y_f.");
disp("Target radius:");
disp(target_radius);
else
    error("Choice of target not understood")
end


%% Initial-state norm and target-pairing error coefficient

% Product L2 norm of the initial state.
if length(y_0_base) == 4*Nx
    % Both components are cubic splines.
    y_01 = y_0_base(1:2*Nx);
    y_02 = y_0_base(2*Nx+1:4*Nx);

    norm_y_0 = sqrt( ...
        iv_cubic_hermite_inner_product(y_01,y_01) ...
        + iv_cubic_hermite_inner_product(y_02,y_02));

elseif length(y_0_base) == 2*Nx-2
    % Both components are P1 functions.
    y_01 = y_0_base(1:Nx-1);
    y_02 = y_0_base(Nx:2*Nx-2);

    norm_y_0 = sqrt( ...
        iv_p1_inner_product(y_01,y_01) ...
        + iv_p1_inner_product(y_02,y_02));

else
    % Concatenated sine coefficients for both components.
    norm_y_0 = norm(intval(y_0_base))/intval('sqrt2');
end


% Target norm used to bound the target-pairing discretization error.
% For spline and P1 targets, the pairing is evaluated directly:
% zero disables this error term; it does not mean the target is zero.
if length(y_f_base) == 4*Nx
    % Cubic spline coefficients.
    norm_y_f = 0;
    y_f1_base = y_f_base(1:2*Nx);

elseif length(y_f_base) == 2*Nx-2
    % P1 nodal coefficients.
    norm_y_f = 0;
    y_f1_base = y_f_base(1:Nx-1);

else
    % Concatenated sine coefficients for both components.
    norm_y_f = norm(intval(y_f_base))/intval('sqrt2');
    y_f1_base = y_f_base(1:length(y_f_base)/2);
end

%% Discretization error bounds

[err,bounds] = iv_coupled_heat_discretization_bound(p_opt_base,norm_y_0,norm_y_f,args);

%% Target support and initial-state contribution

if target_kind == 1
    target_and_initial_term = iv_coupled_partial_target_and_initial(-p_opt, -p_opt_base, y_0_base, y_f1_base, args);
elseif target_kind == 2
    target_and_initial_term = iv_coupled_ball_target_and_initial(-p_opt, -p_opt_base, y_0_base, y_f_base, target_radius, args);
else
    error("Choice of target not understood")
end

%% Computation of support function of control constraints

control_support = iv_coupled_box_control_support(p_opt,args);
disp("Verified control support:");
disp(control_support);

%% Final computations

J_opt = control_support + target_and_initial_term;

bounds.discrete_enclosure = J_opt;

% Record interval uncertainty before adding discretization errors.
bounds.roundoff = rad(J_opt);

disp("Discrete enclosure radius: " + string(bounds.roundoff));

J_opt = J_opt + midrad(0,err);
bounds.final_enclosure = J_opt;

target_not_reachable = isfinite(inf(J_opt)) ...
    && isfinite(sup(J_opt)) && sup(J_opt) < 0;
end