function L_T_star_p = iv_control_adjoint_reference(p, args)
% With C = (I - Dt*A)^(-1), returns B'*(C^(Nt-i+1))'*p at time index i.
% The control-space inner product includes the time-step weight Dt.

    Nt = args{2};
    B_t = transpose(args{5});
    euler_powers = args{8};

    [m,~] = size(B_t);
    
    L_T_star_p = intval(zeros(m,Nt));
    
    
    for i = 1:Nt 
        L_T_star_p(:,i) = B_t*transpose(euler_powers(:,:,Nt-i+2))*p;
    end
end

