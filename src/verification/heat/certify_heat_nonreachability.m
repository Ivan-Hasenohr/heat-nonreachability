function [target_not_reachable,J_opt,control_support,pairing,err,bounds] = certify_heat_nonreachability(p_opt,y_0,y_f_base,args)
% Certify non-reachability for the scalar heat equation with zero initial state.
% This entry point uses positive pointwise control constraints.
% p_opt contains fixed interior nodal values; it need not be an optimizer.
% y_f_base contains cubic Hermite, P1 or sine-series target coefficients.
% The target representation is inferred from its length.
%
% Certification succeeds only when the final enclosure is finite
% and has a strictly negative upper bound.
% bounds records discretization bounds and discrete interval uncertainty.


Dx = args{10};
Nx = args{11};


% Certification requires the verified initializer output.
assert(isa(args{1},'intval') ...
    && isa(args{3},'intval') ...
    && isa(args{4},'intval') ...
    && isa(args{10},'intval'), ...
    'Use verified args from init_heat_fem, not args_float.');


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


assert(isequal(size(args{4}),[Nx-1,Nx-1]), ...
    'The verified evolution matrix has an inconsistent size.');


% The fixed witness contains the interior nodal values.
validateattributes(p_opt,{'double'}, ...
    {'real','finite','column','numel',Nx-1});


%% Select the control constraints

% 1 = positive pointwise bounds; 2 = symmetric L2 bounds.
type_constraints = 1;




% Construct the continuous spline witness from the fixed nodal values.
p_opt_base = iv_natural_cubic_from_nodes(p_opt);

%% Target pairing and coefficient for its discretization error

% A zero norm_y_f disables the target-pairing discretization term
% when the pairing is evaluated directly on the spline witness.
% It does not mean that the target is zero.

initial_data = intval(y_0);
if all(inf(initial_data(:)) == 0) ...
        && all(sup(initial_data(:)) == 0)
    
    if length(y_f_base) == 2*Nx %If y_f is the basis for a cubic spline
        disp("Target is a cubic spline")
        pairing = iv_cubic_hermite_inner_product(y_f_base,p_opt_base); % No discretisation errors on this scalar product <y_f,p_f>


        norm_y_f = 0;
    elseif length(y_f_base) == Nx-1 % Target represented by P1 nodal values.
        pairing = iv_p1_cubic_hermite_inner_product(y_f_base,p_opt_base); % No discretisation errors on this scalar product <y_f,p_f>
        
        norm_y_f = 0; %If y_f is a cubic spline, we have a closed formula for <y_f,p_f> and thus have no discretisation errors! 


    else % Target represented by a sum of sine
        % Include the two zero Dirichlet boundary values.
        pairing =iv_p1_sine_series_pairing([intval(0); intval(p_opt(:)); intval(0)], y_f_base, Dx);

        % L2 norm of the sine-series target, used to bound the target-pairing error.
        norm_y_f = norm(intval(y_f_base))/intval('sqrt2');
    end
else
    error("Only the zero initial state is supported.");

end



%% Compute discretization errors

[err,bounds] = iv_heat_discretization_bound(p_opt_base,y_0,norm_y_f,args,type_constraints);

%% Computing support function:

if type_constraints == 1
    control_support = iv_heat_positive_control_support(p_opt,args);
elseif type_constraints == 2
    control_support = iv_heat_l2_control_support(p_opt,args);
else
    error("Type in a correct constraint index")
end

disp("Verified control support:");
disp(control_support);


%% Finalising the proof:

J_opt = control_support - pairing;
bounds.discrete_enclosure = J_opt;

% Record interval uncertainty before adding discretization errors.
bounds.roundoff = rad(J_opt);

disp("Discrete enclosure radius: " + string(bounds.roundoff));

J_opt = J_opt + midrad(0,err);
bounds.final_enclosure = J_opt;

target_not_reachable = isfinite(inf(J_opt)) ...
    && isfinite(sup(J_opt)) && sup(J_opt) < 0;
end