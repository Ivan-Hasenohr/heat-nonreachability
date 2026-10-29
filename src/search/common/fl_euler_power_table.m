function euler_powers = fl_euler_power_table(Nt, Dt, A)
%FL_EULER_POWER_TABLE Build floating implicit-Euler propagation powers.
% Return an n-by-n-by-(Nt+1) array.
% Slice j+1 approximates (I - Dt*A)^(-j), for j = 0,...,Nt.
% The first slice is the identity.
% These are implicit-Euler powers, not exact matrix exponentials.
% No verified error enclosure is provided.

[~,n] = size(A);
euler_powers = zeros(n,n,Nt+1);
euler_powers(:,:,1) = eye(n);
C = inv(eye(n) - Dt*A);
for i = 2:Nt+1
    euler_powers(:,:,i) = euler_powers(:,:,i-1)*C;
end
end

