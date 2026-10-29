function norm_d1 = iv_cubic_hermite_d1_norm(p_base)

[n,~] = size(p_base);
Nx = fix(n/2);

intval_p_base = intval(p_base);

d1_gram_matrix = zeros(2*Nx,2*Nx);
d1_gram_matrix(1:Nx-1,1:Nx-1) = 36*eye(Nx-1) + diag(-18*ones(Nx-2,1),1) + diag(-18*ones(Nx-2,1),-1);
d1_gram_matrix(Nx:2*Nx,Nx:2*Nx) = 4*eye(Nx+1) + diag(-(1/2)*ones(Nx,1),1) + diag(-(1/2)*ones(Nx,1),-1);
d1_gram_matrix(Nx,Nx) = 2;d1_gram_matrix(2*Nx,2*Nx) = 2;

d1_gram_matrix(1:Nx-1,Nx:2*Nx-2) = diag(-(3/2)*ones(Nx-1,1));d1_gram_matrix(1:Nx-1,Nx+2:2*Nx) = d1_gram_matrix(1:Nx-1,Nx+2:2*Nx) + diag((3/2)*ones(Nx-1,1));
d1_gram_matrix(Nx:2*Nx,1:Nx-1) = transpose(d1_gram_matrix(1:Nx-1,Nx:2*Nx));

d1_gram_matrix = midrad(d1_gram_matrix,0);
d1_gram_matrix = d1_gram_matrix*Nx/15;

norm_d1 = sqrt(intval_p_base'*d1_gram_matrix*intval_p_base);
end