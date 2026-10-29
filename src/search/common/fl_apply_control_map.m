function L_T_u = fl_apply_control_map(u, args_float)
%FL_APPLY_CONTROL_MAP Evaluate the floating discrete control-to-state map.
% u is an m-by-Nt array; B is n-by-m.
% args_float contains floating data from the initializer.
%
% With R = (I - Dt*A)^(-1), return the n-by-1 vector
% L_T_u = Dt * sum_i R^(Nt-i+1) * B * u(:,i).
% Table slice Nt-i+2 represents power Nt-i+1.
% No verified rounding or discretization error enclosure is provided.

    Nt = args_float{2};
    Dt = args_float{3};
    B = args_float{5};
    euler_powers = args_float{8};
    
    [m,~] = size(B);
    L_T_u = zeros(m,1);
    
    for i = 1:Nt
        L_T_u = L_T_u + euler_powers(:,:,Nt-i+2)*B*u(:,i);
    end
    L_T_u = Dt*L_T_u;
end

