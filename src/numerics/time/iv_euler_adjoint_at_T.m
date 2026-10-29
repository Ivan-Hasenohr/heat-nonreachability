function S_T_star_p_f = iv_euler_adjoint_at_T(p_f, args)
% Enclose (I - Dt*A')^(-Nt) * p_f using verified binary powering.
% This is the discrete adjoint evolution at T = Nt*Dt.
% The error relative to the continuous evolution is handled separately.

Nt = args{2};
Dt = args{3};
A = args{4};

[~,n] = size(A);

Ck = verifyinv(eye(n) - Dt*A');

r = floor(log2(Nt)) + 1;

S_T_star_p_f = p_f;

for k = 1:r
    mask = bitget(Nt, k); % Select the power with exponent 2^(k-1).
    if mask
        S_T_star_p_f = Ck * S_T_star_p_f;
    end
    Ck = Ck^2;
end
end