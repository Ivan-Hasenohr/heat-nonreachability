function L_T_u = iv_control_map_reference(u,args)
% Reference implicit-Euler control map.
% With C=(I-Dt*A)^(-1), return Dt*sum_i C^(Nt-i+1)*B*u(:,i).
% args{8} must contain verified powers, including identity in slice 1.
% Continuous-model discretization errors are not included.
    Nt = args{2};
    Dt = args{3};
    B = args{5};
    euler_powers = args{8};

    [m,~] = size(B);
    
    L_T_u = intval(zeros(m,1));
    
    for i = 1:Nt
        L_T_u = L_T_u + euler_powers(:,:,Nt-i+2)*B*intval(u(:,i));
    end

    L_T_u = Dt*L_T_u;
end

