function out = fl_prox_partial_target(p,sigma,y_tilde_T) 
% Y_f = {y_tilde_T} x L^2(0,1)
Nx = fix(length(p)/2)+1;
out = zeros(size(p));

out(1:Nx-1) = p(1:Nx-1) - sigma*y_tilde_T(1:Nx-1);
end