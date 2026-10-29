function pairing = iv_p1_cubic_hermite_inner_product(P1,spline_base)
%Computes the scalar product between a P1 finite element and a cubic spline
%Separates the calculation in two: 
%   - scalar product of P1 and sub-basis of spline defining its values on
%       the mesh
%   - scalar product of P1 and sub-basis of spline defining the values of
%       its derivatives on the mesh


Nx = length(P1)+1;

% "mass" matrix between P1 and spline value basis

toeplitz_val = intval(zeros(Nx-1,1)); toeplitz_val(1) = 14; toeplitz_val(2) = 3;
mat_val = toeplitz(toeplitz_val)/intval('20');

% "mass" matrix between P1 and spline derivative values basis

mat_diff = intval(zeros(Nx-1,Nx+1));
mat_diff(1:Nx-1,1:Nx-1) = eye(Nx-1);
mat_diff(1:Nx-1,3:Nx+1) = mat_diff(1:Nx-1,3:Nx+1) - eye(Nx-1);
mat_diff = mat_diff/intval('30');

% Final computation:

pairing = transpose(P1)*(mat_val*spline_base(1:Nx-1,1) + mat_diff*spline_base(Nx:2*Nx))/intval(Nx);
end