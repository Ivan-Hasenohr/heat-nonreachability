function euler_powers = iv_euler_power_table(Nt, Dt, A)
%IV_EULER_POWER_TABLE Enclose all implicit-Euler propagation powers.
% Inputs: Nt, Dt and A, as defined by the initializer.
% Output: an n-by-n-by-(Nt+1) interval array satisfying
% euler_powers(:,:,j+1) encloses (I - Dt*A)^(-j), for j = 0,...,Nt.
% These are discrete propagation powers, not exact matrix exponentials.
% Storing the full table can require substantial memory.

[~,n] = size(A);
euler_powers = intval(repmat(eye(n),1,1,Nt+1));

Ck = verifyinv(eye(n) - Dt*A);

r = floor(log2(Nt)) + 1;

tic
for k = 1:r
    mask = bitget(0:Nt, k); % Select exponents whose k-th binary digit is one.
    idx = find(mask);
    for j = 1:length(idx)
        euler_powers(:, :, idx(j)) = Ck * euler_powers(:, :,idx(j));
    end
    Ck = Ck^2;
    disp(100*k/r + "%")
end
toc
end

