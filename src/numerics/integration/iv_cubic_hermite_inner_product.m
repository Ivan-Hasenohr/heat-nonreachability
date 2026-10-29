function pairing = iv_cubic_hermite_inner_product(base_1,base_2)

[~,n_spline_1] = size(base_1);
[n,n_spline_2] = size(base_2);
Nx = fix(n/2);

%Basis of the first spline, taking into account 0 values at endpoints

full_coeffs_1 = intval(zeros(2*Nx+2,n_spline_1));
full_coeffs_1(2:Nx,:) = base_1(1:Nx-1,:);
full_coeffs_1(Nx+2:2*Nx+2,:) = base_1(Nx:2*Nx,:);

%Basis of the second spline, taking into account 0 values at endpoints

full_coeffs_2 = intval(zeros(2*Nx+2,n_spline_2));
full_coeffs_2(2:Nx,:) = base_2(1:Nx-1,:);
full_coeffs_2(Nx+2:2*Nx+2,:) = base_2(Nx:2*Nx,:);

%Filling the mass matrix

mass_matrix = zeros(2*(Nx+1),2*(Nx+1));
mass_matrix(1:Nx+1,1:Nx+1) = mass_matrix(1:Nx+1,1:Nx+1) + 78*eye(Nx+1);mass_matrix(1,1) = 39;mass_matrix(Nx+1,Nx+1) = 39;mass_matrix(1:Nx+1,1:Nx+1) = mass_matrix(1:Nx+1,1:Nx+1) + diag(13.5*ones(Nx,1),1) + diag(13.5*ones(Nx,1),-1);
mass_matrix(Nx+2:2*Nx+2,Nx+2:2*Nx+2) = mass_matrix(Nx+2:2*Nx+2,Nx+2:2*Nx+2) + 2*eye(Nx+1);mass_matrix(Nx+2,Nx+2) = 1;mass_matrix(2*Nx+2,2*Nx+2) = 1;mass_matrix(Nx+2:2*Nx+2,Nx+2:2*Nx+2) = mass_matrix(Nx+2:2*Nx+2,Nx+2:2*Nx+2) + diag(-0.75*ones(Nx,1),1) + diag(-0.75*ones(Nx,1),-1);

mass_matrix(1,Nx+2) = 11/2; mass_matrix(Nx+1,2*Nx+2) = -11/2;mass_matrix(1:Nx+1,Nx+2:2*Nx+2) = mass_matrix(1:Nx+1,Nx+2:2*Nx+2) + diag(-(13/4)*ones(Nx,1),1) + diag((13/4)*ones(Nx,1),-1);
mass_matrix(Nx+2:2*Nx+2,1:Nx+1) = transpose(mass_matrix(1:Nx+1,Nx+2:2*Nx+2));

mass_matrix = midrad(mass_matrix,0);
mass_matrix = mass_matrix/(105*Nx);

%Computing the scalar product
pairing = transpose(full_coeffs_1)*mass_matrix*full_coeffs_2;
end