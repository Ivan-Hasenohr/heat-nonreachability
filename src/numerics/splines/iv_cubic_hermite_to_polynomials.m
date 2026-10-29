function polynomial_coeffs = iv_cubic_hermite_to_polynomials(hermite_coeffs)

[n,num_splines] = size(hermite_coeffs);

Nx = fix(n/2);
polynomial_coeffs = intval(zeros(4*Nx,num_splines));
nodal_values = intval(zeros(Nx+1,num_splines));
nodal_values(2:Nx,:) = hermite_coeffs(1:Nx-1,:);
scaled_derivatives = hermite_coeffs(Nx:2*Nx,:);

cell_indices = intval(transpose(0:Nx-1));

polynomial_coeffs(1:4:4*Nx,:) = Nx^3*(2*nodal_values(1:Nx,:)+scaled_derivatives(1:Nx,:) - 2*nodal_values(2:Nx+1,:) + scaled_derivatives(2:Nx+1,:));
polynomial_coeffs(2:4:4*Nx,:) = Nx^2*(-(6*cell_indices+3).*nodal_values(1:Nx,:) - (3*cell_indices+2).*scaled_derivatives(1:Nx,:) + (6*cell_indices+3).*nodal_values(2:Nx+1,:) - (3*cell_indices+1).*scaled_derivatives(2:Nx+1,:));
polynomial_coeffs(3:4:4*Nx,:) = Nx*(6*cell_indices.*(cell_indices+1).*nodal_values(1:Nx,:) + (3*cell_indices+1).*(cell_indices+1).*scaled_derivatives(1:Nx,:) - 6*cell_indices.*(cell_indices+1).*nodal_values(2:Nx+1,:) + cell_indices.*(3*cell_indices+2).*scaled_derivatives(2:Nx+1,:));
polynomial_coeffs(4:4:4*Nx,:) = -(cell_indices+1).^2.*(2*cell_indices-1).*nodal_values(1:Nx,:) - cell_indices.*(cell_indices+1).^2.*scaled_derivatives(1:Nx,:) + cell_indices.^2.*(2*cell_indices+3).*nodal_values(2:Nx+1,:) - cell_indices.^2.*(cell_indices+1).*scaled_derivatives(2:Nx+1,:);

end

