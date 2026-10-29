function target_and_initial_term = iv_coupled_ball_target_and_initial(p_f, p_f_base,y_0_base,y_f_base, target_radius, args)
% Evaluate sigma_Y_f(p_f) - <y_0, S_h_dt^* p_f>.
% Y_f is the product-L2 ball centered at y_f with radius target_radius.
% Its support is <y_f,p_f> + target_radius * ||p_f||.
% S_h_dt^* denotes the discrete adjoint propagator at the final time.
% The caller supplies p_f = -p_opt.

Nx = fix(length(p_f)/2)+1;
target_radius = intval(target_radius);

%% Computation of ||p_f||

% Product L2 norm of the continuous spline witness.
norm_p_f = sqrt( ...
    iv_cubic_hermite_inner_product( ...
        p_f_base(1:2*Nx),p_f_base(1:2*Nx)) ...
    + iv_cubic_hermite_inner_product( ...
        p_f_base(2*Nx+1:4*Nx),p_f_base(2*Nx+1:4*Nx)));
target_and_initial_term = target_radius*norm_p_f;

%% Computation of <y_f, p_f>

if length(y_f_base) == 4*Nx %Both target components are cubic splines
    target_and_initial_term = target_and_initial_term + iv_cubic_hermite_inner_product(y_f_base(1:2*Nx),p_f_base(1:2*Nx)) + iv_cubic_hermite_inner_product(y_f_base(2*Nx+1:end),p_f_base(2*Nx+1:end));
elseif length(y_f_base) == 2*Nx-2 %Both target components are P1 functions.
    target_and_initial_term = target_and_initial_term + iv_p1_cubic_hermite_inner_product(y_f_base(1:Nx-1),p_f_base(1:2*Nx)) + iv_p1_cubic_hermite_inner_product(y_f_base(Nx:end),p_f_base(2*Nx+1:end));
else %Both targets are sums of sinuses (of same length!)
    n = length(y_f_base); 
    if rem(n,2)
        error("Both target components must have the same number of sine coefficients.")
    end
    num_modes = n/2;

    
    target_and_initial_term = target_and_initial_term + iv_p1_sine_series_pairing_vectorized([0;p_f(1:Nx-1);0], y_f_base(1:num_modes), 1/intval(Nx)) + iv_p1_sine_series_pairing_vectorized([0;p_f(Nx:end);0], y_f_base(num_modes+1:end), 1/intval(Nx));
end

%% Initial-state pairing with the discrete adjoint evolution

% Apply the discrete adjoint propagator to the witness.
adjoint_at_T = iv_euler_adjoint_at_T(p_f, args);


if length(y_0_base) == 4*Nx % If y_0 is two cubic splines
    target_and_initial_term = target_and_initial_term - iv_p1_cubic_hermite_inner_product(adjoint_at_T(1:Nx-1),y_0_base(1:2*Nx)) - iv_p1_cubic_hermite_inner_product(adjoint_at_T(Nx:2*Nx-2),y_0_base(2*Nx+1:4*Nx));

elseif length(y_0_base) == 2*Nx-2 % If y_0 is two finite elements
    target_and_initial_term = target_and_initial_term - iv_p1_inner_product(adjoint_at_T(1:Nx-1),y_0_base(1:Nx-1)) - iv_p1_inner_product(adjoint_at_T(Nx:2*Nx-2),y_0_base(Nx:2*Nx-2));

else
    n = length(y_0_base); 
    if rem(n,2)
        error("If y_0 is a sum of sin, then it must have the same length of decomposition for y_01 and y_02.")
    end
    num_modes = n/2;
    target_and_initial_term = target_and_initial_term - iv_p1_sine_series_pairing_vectorized([0;adjoint_at_T(1:Nx-1);0], y_0_base(1:num_modes), 1/intval(Nx))  - iv_p1_sine_series_pairing_vectorized([0;adjoint_at_T(Nx:2*Nx-2);0], y_0_base(num_modes+1:end), 1/intval(Nx));

end

end