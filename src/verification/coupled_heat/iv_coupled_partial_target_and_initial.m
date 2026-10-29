function target_and_initial_term = iv_coupled_partial_target_and_initial(p_f, p_f_base,y_0_base,y_f1_base, args)
% Evaluate sigma_Y_f(p_f) - <y_0, S_h_dt^* p_f>.
% Y_f = {y_f1} x L2(0,1); its support is infinite unless p_f2 is zero.
% S_h_dt^* denotes the discrete adjoint propagator at the final time.
% The caller supplies p_f = -p_opt.

Nx = fix(length(p_f)/2)+1;

% The unconstrained target component requires an exactly zero witness component.
if any(p_f(Nx:2*Nx-2) ~= 0)
    target_and_initial_term = Inf;
else
    %% Computation of <y_f, p_f>

    if length(y_f1_base) == 2*Nx %y_f1 is a cubic spline
        target_and_initial_term = iv_cubic_hermite_inner_product(y_f1_base,p_f_base(1:2*Nx));
    
    elseif length(y_f1_base) == Nx-1 %y_f1 is a P1 finite element
        target_and_initial_term = iv_p1_cubic_hermite_inner_product(y_f1_base,p_f_base(1:2*Nx));
    
    else %y_f1 is a sum of sin
        % Integrate the first component, including its zero boundary values.
        target_and_initial_term = iv_p1_sine_series_pairing([intval(0); intval(p_f(1:Nx-1)); intval(0)], y_f1_base, 1/intval(Nx));
    end
    
%% Initial-state pairing with the discrete adjoint evolution

    adjoint_at_T = iv_euler_adjoint_at_T(p_f, args);

    % Subtract the initial-state pairing because the caller passes p_f = -p_opt.
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
        target_and_initial_term = target_and_initial_term - iv_p1_sine_series_pairing([intval(0); intval(adjoint_at_T(1:Nx-1)); intval(0)], y_0_base(1:num_modes), 1/intval(Nx));
        
        target_and_initial_term = target_and_initial_term - iv_p1_sine_series_pairing([intval(0); intval(adjoint_at_T(Nx:2*Nx-2)); intval(0)], y_0_base(num_modes+1:end), 1/intval(Nx));

    end
end
end
