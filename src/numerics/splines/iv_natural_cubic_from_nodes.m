function hermite_coeffs = iv_natural_cubic_from_nodes(nodal_values)
    Nx = length(nodal_values) + 1;

    hermite_coeffs = intval(zeros(2*Nx,1));
    hermite_coeffs(1:Nx-1) = nodal_values;

    coupling_first_row = intval(zeros(Nx-1,1)); coupling_first_row(2) =  -6;
    coupling_first_column = intval(zeros(Nx-1,1)); coupling_first_column(2) =  6;
    nodal_values_d1_coupling(2:Nx,:) = toeplitz(coupling_first_column,coupling_first_row);
    nodal_values_d1_coupling(1,1) = -6; nodal_values_d1_coupling(Nx+1,Nx-1) = 6; 
    
    d1_system_first_column = intval(zeros(Nx+1,1)); 
    d1_system_first_column(1) = 8; d1_system_first_column(2) = 2;
    d1_system_matrix = toeplitz(d1_system_first_column);
    d1_system_matrix(1,1) = 4; d1_system_matrix(Nx+1,Nx+1) = 4;
    
    hermite_coeffs(Nx:2*Nx) = verifylss(d1_system_matrix,-nodal_values_d1_coupling*nodal_values);

    
