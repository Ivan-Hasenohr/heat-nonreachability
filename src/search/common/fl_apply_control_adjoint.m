function L_T_star_p = fl_apply_control_adjoint(p, args_float)
%FL_APPLY_CONTROL_ADJOINT Evaluate the floating discrete control adjoint.
% p is an n-by-1 vector; B is n-by-m.
% Return an m-by-Nt array whose column i is
% B' * (R')^(Nt-i+1) * p, with R = (I - Dt*A)^(-1).
% Table slice Nt-i+2 represents power Nt-i+1.
%
% The control-space pairing includes the time weight Dt:
% p'*L_T(u) = Dt * sum_i u(:,i)' * L_T_star(p)(:,i),
% up to floating-point errors.
% No verified rounding or discretization error enclosure is provided.
    Nt = args_float{2};
    B_t = transpose(args_float{5});
    euler_powers = args_float{8};

    [m,~] = size(B_t);
    L_T_star_p = zeros(m,Nt);
    for i = 1:Nt 
        L_T_star_p(:,i) = B_t*transpose(euler_powers(:,:,Nt-i+2))*p;
    end

end

