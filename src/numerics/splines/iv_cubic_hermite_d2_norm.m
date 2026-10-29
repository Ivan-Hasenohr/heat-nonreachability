function norm_d2 = iv_cubic_hermite_d2_norm(p_base)

[n,n_spline] = size(p_base);
Nx = fix(n/2);

d2_gram_matrix = zeros(2*(Nx+1),2*(Nx+1));
d2_gram_matrix(1:Nx+1,1:Nx+1) = 24*eye(Nx+1) + diag(-12*ones(Nx,1),1) + diag(-12*ones(Nx,1),-1); d2_gram_matrix(1,1) = 12;d2_gram_matrix(Nx+1,Nx+1) = 12;
d2_gram_matrix(Nx+2:2*Nx+2,Nx+2:2*Nx+2) = 8*eye(Nx+1)+ diag(2*ones(Nx,1),1) + diag(2*ones(Nx,1),-1);d2_gram_matrix(Nx+2,Nx+2) = 4;d2_gram_matrix(2*Nx+2,2*Nx+2) = 4;

d2_gram_matrix(1,Nx+2) = 6; d2_gram_matrix(Nx+1,2*Nx+2) = -6;d2_gram_matrix(1:Nx+1,Nx+2:2*Nx+2) = d2_gram_matrix(1:Nx+1,Nx+2:2*Nx+2) + diag(6*ones(Nx,1),1) + diag(-6*ones(Nx,1),-1);
d2_gram_matrix(Nx+2:2*Nx+2,1:Nx+1) = transpose(d2_gram_matrix(1:Nx+1,Nx+2:2*Nx+2));


d2_gram_matrix = midrad(d2_gram_matrix,0);
d2_gram_matrix = d2_gram_matrix*Nx^3;

%Adapt basis to full basis (including zero boundary conditions)
full_coeffs = intval(zeros(2*Nx+2,n_spline));
full_coeffs(2:Nx,:) = p_base(1:Nx-1,:);
full_coeffs(Nx+2:2*Nx+2,:) = p_base(Nx:2*Nx,:);


norm_d2 = sqrt(transpose(full_coeffs)*d2_gram_matrix*full_coeffs);
end