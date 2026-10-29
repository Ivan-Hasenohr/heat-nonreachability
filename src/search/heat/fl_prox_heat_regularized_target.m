function out = fl_prox_heat_regularized_target(p, sigma,y_tilde_T, lambda_reg,args_float)
% Proximal map of <y_tilde_T,p> + (lambda_reg/2)*||A'*p||_2^2.
A = args_float{4};
[n,~] = size(A);
out = (eye(n) + sigma*lambda_reg*A*transpose(A))\(p - sigma*y_tilde_T);
end

