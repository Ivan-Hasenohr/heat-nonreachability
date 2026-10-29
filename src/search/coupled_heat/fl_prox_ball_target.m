function out = fl_prox_ball_target(p,sigma,y_tilde_T, target_radius) 
% Y_f = {y, \|y-y_tilde_T\| <= target_radius}

Nx = fix(length(p)/2)+1;

shifted_p = p - sigma*y_tilde_T;
shifted_p_norm = norm(shifted_p)/sqrt(Nx);

out = max(0,1-sigma*target_radius/shifted_p_norm)*shifted_p;
end