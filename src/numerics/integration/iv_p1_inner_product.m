function pairing = iv_p1_inner_product(p,q)
Nx = length(p) + 1;

intval_int_range = intval(linspace(1,Nx-1,Nx-1));

intval_D_diag =( 1 + 0.5*cos(intval('pi')/intval(Nx)*intval_int_range))/((3/2)*intval(Nx));

intval_V = sqrt(2/intval(Nx))*sin(intval('pi')/intval(Nx)*intval_int_range.*transpose(intval_int_range));

intval_M = intval_V*diag(intval_D_diag)*intval_V;

pairing = p'*intval_M*q;
end